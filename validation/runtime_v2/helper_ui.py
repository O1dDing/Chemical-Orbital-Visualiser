# SPDX-License-Identifier: CPAL-1.0
#
# The contents of this file are subject to the Common Public Attribution License
# Version 1.0 (the “License”); you may not use this file except in compliance with the
# License. You may obtain a copy of the License at
# https://opensource.org/license/CPAL-1.0. The License is based on the Mozilla Public
# License Version 1.1 but Sections 14 and 15 have been added to cover use of software
# over a computer network and provide for limited attribution for the Original
# Developer. In addition, Exhibit A has been modified to be consistent with Exhibit B.
#
# Software distributed under the License is distributed on an “AS IS” basis, WITHOUT
# WARRANTY OF ANY KIND, either express or implied. See the License for the specific
# language governing rights and limitations under the License.
#
# The Original Code is Chemical Orbital Visualiser (COV).
#
# The Original Developer is O1dDing, who is also the Initial Developer.
#
# The Initial Developer of the Original Code is O1dDing. All portions of the code
# written by O1dDing are Copyright (c) 2026 O1dDing. All Rights Reserved.
#
# Contributor O1dDing.

"""Visible COV Helper. Startup is read-only and never imports or resumes work."""
from pathlib import Path
import argparse
import json
import os
import queue
import threading
import tkinter as tk
from tkinter import messagebox, ttk

from helper_backend import HELPER_VERSION, HelperBackend


class HelperWindow:
    def __init__(self, root, backend):
        self.root, self.backend = root, backend
        self.busy = False
        self.results = queue.Queue()
        self.closed = False
        self.root.title('COV Helper ' + HELPER_VERSION)
        self.root.geometry('1040x740')
        self.root.minsize(880, 650)
        style = ttk.Style(root)
        style.theme_use('clam')
        style.configure('.', font=('Microsoft YaHei UI', 10))
        style.configure('Title.TLabel', font=('Microsoft YaHei UI', 22, 'bold'))
        style.configure('Count.TLabel', font=('Segoe UI', 25, 'bold'))
        style.configure('Quiet.TLabel', foreground='#536274')
        main = ttk.Frame(root, padding=24)
        main.pack(fill='both', expand=True)
        main.columnconfigure(0, weight=1)
        main.rowconfigure(5, weight=1)
        ttk.Label(main, text='COV Helper  ' + HELPER_VERSION, style='Title.TLabel').grid(row=0, column=0, sticky='w')
        ttk.Label(main, text='自适应物理核运行器 3.0  ·  REF-001 计算管理', style='Quiet.TLabel').grid(row=1, column=0, sticky='w', pady=(2, 16))
        self.banner = tk.StringVar(value='正在读取状态；计算控制未启用')
        ttk.Label(main, textvariable=self.banner, font=('Microsoft YaHei UI', 12, 'bold')).grid(row=2, column=0, sticky='w')
        cards = ttk.Frame(main)
        cards.grid(row=3, column=0, sticky='ew', pady=16)
        self.counts = {}
        for column, (key, label) in enumerate((('candidate_collected', '候选收集完成'), ('partial', '部分完成'), ('needs_review', '待复核'), ('not_started', '未开始'))):
            cards.columnconfigure(column, weight=1)
            frame = ttk.LabelFrame(cards, text=label, padding=(15, 8))
            frame.grid(row=0, column=column, sticky='ew', padx=(0, 12 if column < 3 else 0))
            value = tk.StringVar(value='—')
            self.counts[key] = value
            ttk.Label(frame, textvariable=value, style='Count.TLabel').pack(anchor='w')
        toolbar = ttk.Frame(main)
        toolbar.grid(row=4, column=0, sticky='ew', pady=(0, 12))
        self.refresh_button = ttk.Button(toolbar, text='刷新状态', command=self.refresh)
        self.refresh_button.pack(side='left')
        self.check_button = ttk.Button(toolbar, text='检查部署（不计算）', command=lambda: self.submit('check', backend.check))
        self.check_button.pack(side='left', padx=8)
        ttk.Button(toolbar, text='打开计算目录', command=lambda: os.startfile(backend.work)).pack(side='left')
        self.detail = tk.Text(main, wrap='word', height=10, font=('Microsoft YaHei UI', 10),
                              relief='solid', borderwidth=1, padx=12, pady=10, background='#ffffff')
        self.detail.grid(row=5, column=0, sticky='nsew')
        controls = ttk.LabelFrame(main, text='手动计算操作', padding=12)
        controls.grid(row=6, column=0, sticky='ew', pady=(16, 0))
        self.enabled = tk.BooleanVar(value=False)
        self.enable_box = ttk.Checkbutton(controls, text='允许本窗口手动控制计算（每次打开默认关闭）', variable=self.enabled, command=self.toggle)
        self.enable_box.grid(row=0, column=0, columnspan=8, sticky='w', pady=(0, 10))
        self.limit, self.workers, self.cases = tk.StringVar(value='1'), tk.StringVar(value='1'), tk.StringVar()
        ttk.Label(controls, text='案例上限').grid(row=1, column=0, sticky='w')
        ttk.Spinbox(controls, from_=1, to=273, textvariable=self.limit, width=5).grid(row=1, column=1, padx=7)
        ttk.Label(controls, text='并行槽位').grid(row=1, column=2)
        ttk.Combobox(controls, values=('1', '2'), textvariable=self.workers, state='readonly', width=4).grid(row=1, column=3, padx=7)
        ttk.Label(controls, text='案例 ID（留空按队列）').grid(row=1, column=4)
        ttk.Entry(controls, textvariable=self.cases, width=25).grid(row=1, column=5, columnspan=3, padx=7, sticky='ew')
        self.action_buttons = []
        actions = [('开始批次', self.start), ('内存暂停', lambda: self.control('hold')), ('恢复计算', lambda: self.control('resume')),
                   ('阶段后暂停', lambda: self.control('pause')), ('停止并冷保存', lambda: self.control('interrupt')),
                   ('阶段后退出', lambda: self.control('shutdown'))]
        for column, (label, callback) in enumerate(actions):
            button = ttk.Button(controls, text=label, command=callback, state='disabled')
            button.grid(row=2, column=column, padx=(0, 6), pady=(12, 0))
            self.action_buttons.append(button)
        ttk.Label(main, text='打开或关闭此面板不会启动、恢复或停止计算。候选收集完成不等于科学验收通过。', style='Quiet.TLabel', wraplength=970).grid(row=7, column=0, sticky='w', pady=(12, 0))
        self.root.protocol('WM_DELETE_WINDOW', self.close)
        self.root.after(50, self.poll)
        self.refresh()

    def toggle(self):
        self.backend.arm(self.enabled.get())
        self.sync_buttons()

    def sync_buttons(self):
        for button in self.action_buttons:
            button.configure(state='normal' if self.backend.armed and not self.busy else 'disabled')
        self.refresh_button.configure(state='disabled' if self.busy else 'normal')
        self.check_button.configure(state='disabled' if self.busy else 'normal')

    def submit(self, kind, operation):
        if self.busy or self.closed:
            return
        self.busy = True
        self.sync_buttons()
        def work():
            try:
                self.results.put((kind, operation(), None))
            except Exception as error:
                self.results.put((kind, None, str(error)))
        threading.Thread(target=work, daemon=True).start()

    def refresh(self):
        self.submit('status', self.backend.status)

    def show_text(self, text):
        self.detail.configure(state='normal')
        self.detail.delete('1.0', 'end')
        self.detail.insert('1.0', text)
        self.detail.configure(state='disabled')

    def poll(self):
        if self.closed:
            return
        try:
            kind, result, error = self.results.get_nowait()
        except queue.Empty:
            self.root.after(100, self.poll)
            return
        self.busy = False
        self.sync_buttons()
        if error:
            self.banner.set('操作未完成；请查看下方信息')
            self.show_text(error)
        elif kind == 'status':
            running = bool(result.get('work_lease_held'))
            mode = (result.get('control') or {}).get('mode', 'unknown')
            translated = {'pause': '阶段暂停', 'hold': '内存暂停', 'run': '运行', 'interrupt': '停止并保存', 'shutdown': '阶段后退出'}.get(mode, mode)
            self.banner.set(('协调器运行中' if running else '空闲 · 无运行中的协调器') + '    |    控制状态：' + translated)
            for key, value in self.counts.items():
                value.set(str(result.get('counts', {}).get(key, 0)))
            binding = result.get('binding', {})
            text = ('启动方式：只读取状态；未自动提交、导入、启动或恢复计算。\n'
                    f'运行器身份：{binding.get("runtime_identity", "未绑定")}\n'
                    f'计算目录：{self.backend.work}\n'
                    f'共享配置：{self.backend.config_path}\n'
                    '已有计算与恢复记录保持原位；候选状态不表示科学验收结果。\n'
                    '内存暂停保留 RAM 现场；需要关机时须先取得冷保存回执。')
            text += ('\n资源策略：16物理核总预算14，单任务4/5/7核；'
                     '8物理核总预算7，单任务2/3/4核。最多两个任务，内存暂停继续占用原配额。')
            if result.get('progress_text'):
                text += '\n\n' + result['progress_text']
            self.show_text(text)
        elif kind == 'check':
            self.banner.set('部署检查通过 · 没有启动计算')
            self.show_text(f'冻结输入：{result["verified_inputs"]} 个\n物理核总预算：{result["core_budget"]}\ncalculation_started: false\n\n' + json.dumps(result, ensure_ascii=False, indent=2))
        else:
            self.show_text(json.dumps(result, ensure_ascii=False, indent=2))
            self.root.after(200, self.refresh)
        self.root.after(100, self.poll)

    def start(self):
        values = (self.limit.get(), self.workers.get(), self.cases.get())
        self.submit('start', lambda: self.backend.start(*values))

    def control(self, command):
        if command == 'interrupt' and not messagebox.askokcancel('停止并冷保存', '结束所属计算并尝试保存 CHK/RWF；未保存的阶段工作可能丢失。继续？', parent=self.root):
            return
        self.submit('control', lambda: self.backend.control(command))

    def close(self):
        self.closed = True
        self.backend.arm(False)
        self.root.destroy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--home', type=Path, required=True)
    parser.add_argument('--check', action='store_true', help='Read-only deployment check without opening a window')
    args = parser.parse_args()
    backend = HelperBackend(args.home)
    if args.check:
        print(json.dumps({'helper_version': HELPER_VERSION, 'manual_controls_enabled': backend.armed,
                          'status': backend.status(), 'deployment': backend.check()}, ensure_ascii=False, indent=2))
        return
    from state_store import lease
    with lease(args.home / 'helper-operation.lock'):
        root = tk.Tk()
        HelperWindow(root, backend)
        root.mainloop()


if __name__ == '__main__':
    main()
