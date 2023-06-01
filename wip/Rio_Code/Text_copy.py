def write_excel_xlsx(path, sheet_name, value):
    import numpy as np
    index = len(value)
    workbook = openpyxl.Workbook()  
    sheet = workbook.active 
    sheet.title = sheet_name  
    for i in range(0, index):
        for j in range(0, len(value[i])):
            sheet.cell(row=i + 1, column=j + 1, value=str(value[i][j])) 
    workbook.save(path)  
    print("Successfully create a excel file")