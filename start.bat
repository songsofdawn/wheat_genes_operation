@echo off
title start running streamlit app...
:: 检查是否在当前目录下找到 app.py
if not exist app.py (
    echo [错误] 未能在当前目录下找到 app.py 文件！
    pause
    exit
)

:: 运行 Streamlit
streamlit run app.py

pause