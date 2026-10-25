@echo off
REM SPDX-License-Identifier: CPAL-1.0
REM
REM The contents of this file are subject to the Common Public Attribution License
REM Version 1.0 (the “License”); you may not use this file except in compliance with the
REM License. You may obtain a copy of the License at
REM https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
REM License Version 1.1 but Sections 14 and 15 have been added to cover use of software
REM over a computer network and provide for limited attribution for the Original
REM Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
REM
REM Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
REM WARRANTY OF ANY KIND, either express or implied. See the License for the specific
REM language governing rights and limitations under the License.
REM
REM The Original Code is Chemical Orbital Visualiser (COV).
REM
REM The Original Developer is O1dDing, who is also the Initial Developer.
REM
REM The Initial Developer of the Original Code is O1dDing. All portions of the code
REM written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
REM
REM Contributor O1dDing.

setlocal
chcp 65001 >nul
cd /d "%~dp0"
set "COV_CONFIG=%~dp0config.json"
if exist "%~dp0..\..\..\..\active.json" set "COV_CONFIG=%~dp0..\..\..\..\config.json"
if not exist "%COV_CONFIG%" (
  copy /y "config.example.json" "%COV_CONFIG%" >nul
  echo Created config.json. Check its three paths before choosing Import or Run.
)
if exist "F:\Dev\Python312\python.exe" (
  "F:\Dev\Python312\python.exe" -X utf8 "%~dp0resume.py" --config "%COV_CONFIG%" menu
) else (
  where py >nul 2>nul
  if errorlevel 1 (
    python -X utf8 "%~dp0resume.py" --config "%COV_CONFIG%" menu
  ) else (
    py -3.12 -X utf8 "%~dp0resume.py" --config "%COV_CONFIG%" menu
  )
)
if errorlevel 1 pause
endlocal
