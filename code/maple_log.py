import csv
import openpyxl
import os

def process_maple_log(input_log, output_csv, output_xlsx):
    if not os.path.exists(input_log):
        print(f"❌ 错误：找不到文件 {input_log}")
        return

    # 核心配置
    target_phrase = "The system has given number of real solution(s)"
    all_results = []       # 用于 Excel：[索引, 状态, 完整内容]
    promising_indices = [] # 用于 CSV：[索引]

    current_index = None
    current_block = []

    print(f"正在读取日志并生成文件...")

    with open(input_log, 'r', encoding='utf-8') as f:
        for line in f:
            text = line.strip()
            if not text or text.startswith("---"):
                continue

            # 识别索引行
            if text.isdigit():
                # 处理上一个块
                if current_index is not None:
                    block_text = "\n".join(current_block)
                    is_ok = target_phrase in block_text
                    status = "Promising" if is_ok else "No/Error"
                    
                    all_results.append([current_index, status, block_text])
                    if is_ok:
                        promising_indices.append(current_index)

                # 重置进入新块
                current_index = text
                current_block = []
            else:
                current_block.append(text)

        # 处理最后一个块
        if current_index:
            block_text = "\n".join(current_block)
            is_ok = target_phrase in block_text
            status = "Promising" if is_ok else "No/Error"
            all_results.append([current_index, status, block_text])
            if is_ok:
                promising_indices.append(current_index)

    # --- 输出 1: CSV 文件 ---
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        for idx in promising_indices:
            writer.writerow([idx])

    # --- 输出 2: Excel 文件 ---
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Analysis"
    ws.append(["Index", "Status", "Full Maple Content"]) # 表头
    
    for row in all_results:
        ws.append(row)

    # 简单美化 Excel
    ws.column_dimensions['B'].width = 15
    ws.column_dimensions['C'].width = 80
    wb.save(output_xlsx)

    print(f" 处理完成！")
    print(f" 总计处理：{len(all_results)} 个网络")
    print(f" 发现：{len(promising_indices)} 个")
    print(f" CSV 文件：{output_csv}")
    print(f" Excel 文件：{output_xlsx}")

if __name__ == "__main__":
    # 配置路径
    log_path = "D:/multistability/maple_output3-8_log.txt"
    csv_path = "D:/multistability/candidate_indices.csv"
    xlsx_path = "D:/multistability/maple_results_full_check.xlsx"

    process_maple_log(log_path, csv_path, xlsx_path)