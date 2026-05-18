@echo off
echo ====================================
echo  GitHub 项目初始化脚本
echo ====================================
echo.

REM 检查git是否安装
where git >nul 2>&1
if %errorlevel% neq 0 (
    echo [错误] 未检测到Git，请先安装 Git for Windows
    echo 下载地址: https://git-scm.com/download/win
    echo.
    echo 安装后重新运行此脚本
    pause
    exit /b 1
)

echo [OK] 检测到Git
echo.

REM 初始化仓库
echo [1/4] 初始化Git仓库...
git init
if %errorlevel% neq 0 (
    echo [错误] Git初始化失败
    pause
    exit /b 1
)

REM 添加文件
echo [2/4] 添加文件到暂存区...
git add .

REM 提交
echo [3/4] 提交到本地仓库...
git commit -m "Initial commit - protein linker analysis toolkit"

REM 显示状态
echo [4/4] 显示仓库状态...
git status

echo.
echo ====================================
echo  本地仓库初始化完成！
echo ====================================
echo.
echo 下一步操作:
echo 1. 登录 GitHub: https://github.com
echo 2. 创建新仓库 (New Repository)
echo 3. 将仓库名称设置为: analyze_linker_tool
echo 4. 复制仓库的HTTPS地址
echo 5. 运行以下命令连接远程仓库:
echo.
echo    git remote add origin [你的仓库地址]
echo    git branch -M main
echo    git push -u origin main
echo.
echo 或者使用 GitHub Desktop 克隆此仓库
echo.
pause
