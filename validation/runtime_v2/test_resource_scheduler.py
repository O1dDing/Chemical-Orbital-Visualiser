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

"""Admission races, durable crash reconciliation and full-tier fairness."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tempfile
import threading
import unittest
import os
import subprocess
import sys

from policy import core_budget, preferred_cores, grants
from resource_scheduler import Ledger, AdmissionDenied
from state_store import read_json


class ResourceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / 'ledger.json'
        self.living = {(11, 101)}
        self.jobs = set()

    def ledger(self, physical=16, owner=None):
        return Ledger(self.path, {'physical_cores': physical}, owner or {'pid':11,'creation_filetime':101},
            lambda p: (p['pid'],p['creation_filetime']) in self.living, lambda n: n in self.jobs)

    def test_tiers_and_thresholds(self):
        for n, sixteen, eight in [(1,4,2),(299,4,2),(300,5,3),(449,5,3),(450,7,4),(574,7,4),(5000,7,4)]:
            for shell in (False, True):
                self.assertEqual(preferred_cores(n,shell,16),sixteen)
                self.assertEqual(preferred_cores(n,shell,8),eight)

    def test_all_required_sixteen_pairs_and_third_job(self):
        ledger=self.ledger()
        for a,b in [(4,4),(4,5),(5,5),(7,4),(7,5),(7,7)]:
            x=ledger.reserve(a,label='a'); y=ledger.reserve(b,label='b')
            with self.assertRaises(AdmissionDenied): ledger.reserve(1,label='third')
            y.release();x.release()
        x=ledger.reserve(7,label='heavy')
        with self.assertRaises(AdmissionDenied): ledger.reserve(8,label='over14')
        x.release()

    def test_eight_four_plus_three_but_not_four_plus_four(self):
        ledger=self.ledger(8); x=ledger.reserve(4,label='heavy')
        with self.assertRaises(AdmissionDenied): ledger.reserve(4,label='heavy2')
        y=ledger.reserve(3,label='medium');y.release();x.release()
        self.assertEqual(core_budget(8),7)

    def test_atomic_racing_reservations(self):
        barrier=threading.Barrier(12)
        def contender(i):
            barrier.wait()
            try: return self.ledger().reserve(7,label=str(i))
            except AdmissionDenied: return None
        with ThreadPoolExecutor(max_workers=12) as pool:
            admitted=[x for x in pool.map(contender,range(12)) if x]
        self.assertEqual(len(admitted),2)
        self.assertEqual(sum(r['allocated'] for r in read_json(self.path)['active'].values()),14)
        for x in admitted:x.release()

    def test_ram_pause_retains_and_readmit_after_exit(self):
        ledger=self.ledger(); x=ledger.reserve(7,label='heavy')
        ledger.update_tree(x.token,{'job_name':'job','process':{'pid':22,'creation_filetime':102}})
        ledger.hold(x.token,True)
        self.assertEqual(read_json(self.path)['active'][x.token]['allocated'],7)
        with self.assertRaises(AdmissionDenied):x.release()
        with self.assertRaises(AdmissionDenied):ledger.update_tree(x.token,{'job_name':'concurrent-probe'})
        ledger.hold(x.token,False);ledger.update_tree(x.token,None);x.release()
        y=ledger.reserve(7,label='resume');y.release()

    def test_crash_does_not_forget_living_child_or_job(self):
        ledger=self.ledger();x=ledger.reserve(7,label='old')
        ledger.update_tree(x.token,{'job_name':'orphan','process':{'pid':22,'creation_filetime':102}})
        self.living.clear();self.jobs.add('orphan')
        replacement=self.ledger(owner={'pid':33,'creation_filetime':103})
        self.living.add((33,103))
        y=replacement.reserve(7,label='new')
        with self.assertRaises(AdmissionDenied):replacement.reserve(1,label='third')
        y.release();self.jobs.clear();self.living.add((22,102))
        y=replacement.reserve(7,label='root-still-alive'); y.release()
        self.living.remove((22,102))
        # PID reuse is not the old writer; native creation identity must agree.
        self.living.add((22,999))
        y=replacement.reserve(7,label='reconciled')
        self.assertEqual(len(read_json(self.path)['active']),1)
        self.assertTrue(any(r['event']=='crash_reconciled_release' for r in read_json(self.path)['events']))
        y.release()

    def test_fifo_wait_never_waterfills_or_bypasses_heavy(self):
        self.assertEqual(grants([7,4],6),[0,0])
        self.assertEqual(grants([7,4],14),[7,4])
        self.assertEqual(grants([4,4],7),[4,0])

    def test_cross_process_atomic_race(self):
        gate = Path(self.temp.name) / 'go'
        script = Path(self.temp.name) / 'race.py'
        script.write_text('''import sys,time
from pathlib import Path
sys.path.insert(0,sys.argv[1])
from resource_scheduler import Ledger,AdmissionDenied
while not Path(sys.argv[3]).exists():time.sleep(.01)
ledger=Ledger(sys.argv[2],{'physical_cores':16},{'pid':11,'creation_filetime':101},lambda p:True)
try:
 ledger.reserve(7,label=sys.argv[4]);print('admitted')
except AdmissionDenied:print('denied')
''')
        children = [subprocess.Popen([sys.executable, '-B', str(script), str(Path(__file__).parent),
                    str(self.path), str(gate), str(i)], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0)) for i in range(6)]
        gate.write_text('go')
        try:
            results = [p.communicate(timeout=10) for p in children]
            self.assertEqual([p.returncode for p in children],[0]*6)
            self.assertEqual(sum(out.strip()==b'admitted' for out,err in results),2)
            self.assertEqual(sum(r['allocated'] for r in read_json(self.path)['active'].values()),14)
        finally:
            for p in children:
                if p.poll() is None:p.kill();p.wait()

    @unittest.skipUnless(os.name=='nt','Native crash cleanup requires Windows')
    def test_native_crash_at_creation_reconciles_without_losing_child(self):
        import windows_job
        script=Path(self.temp.name)/'crash.py'
        script.write_text('''import sys,os
from pathlib import Path
sys.path.insert(0,sys.argv[1]);sys.path.append(str(Path(sys.argv[1]).parent/'paused-20260906/runner-source'))
import windows_job
from policy import core_budget
folder=Path(sys.argv[2]);(folder/'scratch').mkdir()
cores=min(7,core_budget(windows_job.host_info()['physical_cores']))
windows_job.run_tree([sys.executable,'-c','import time;time.sleep(60)'],folder,cores,1,60,on_started=lambda value:os._exit(0))
''')
        child=subprocess.run([sys.executable,'-B',str(script),str(Path(__file__).parent),self.temp.name],
            env=dict(os.environ,LOCALAPPDATA=self.temp.name),capture_output=True,timeout=10,
            creationflags=subprocess.CREATE_NO_WINDOW)
        self.assertEqual(child.returncode,0,child.stderr)
        path=Path(self.temp.name)/'COV/adaptive-physical-resources.json'
        old=next(iter(read_json(path)['active'].values()))
        self.assertFalse(windows_job.process_alive(old['owner']))
        self.assertFalse(windows_job.job_alive(old['tree']['job_name']))
        if old['tree']['process']:
            self.assertFalse(windows_job.process_alive(old['tree']['process']))
        ledger=Ledger(path,{'physical_cores':16},windows_job.process_identity(),windows_job.process_alive,windows_job.job_alive)
        x=ledger.reserve(7,label='after-crash');x.release()
        self.assertEqual(read_json(path)['active'],{})


if __name__=='__main__':unittest.main(verbosity=2)
