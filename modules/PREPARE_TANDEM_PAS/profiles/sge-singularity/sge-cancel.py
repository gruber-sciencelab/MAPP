#!/usr/bin/env python3
import subprocess as sp
import shlex
import sys

jobid_list = ', '.join(sys.argv[1:])

sp.check_call(shlex.split(f"qdel {jobid_list}"))

