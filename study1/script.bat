@echo off

set lambda=1
set v0=1
set r0=1
set tw=1.5708
set simDur=2e5
set dT=0.01
set stageOrTime=1
set log_stages=0

REM a = 0.2
START FOR %%a IN (0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99) DO  (..\sim.exe %lambda% %v0% %r0% %tw% %simDur% %dT% a e 0.2 %%a %stageOrTime% %log_stages% 0)

REM a = 0.4
START FOR %%a IN (0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99) DO  (..\sim.exe %lambda% %v0% %r0% %tw% %simDur% %dT% a e 0.4 %%a %stageOrTime% %log_stages% 0)

REM a = 0.6
START FOR %%a IN (0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99) DO  (..\sim.exe %lambda% %v0% %r0% %tw% %simDur% %dT% a e 0.6 %%a %stageOrTime% %log_stages% 0)

REM a = 0.8
START FOR %%a IN (0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99) DO  (..\sim.exe %lambda% %v0% %r0% %tw% %simDur% %dT% a e 0.8 %%a %stageOrTime% %log_stages% 0)

REM a = 1
START FOR %%a IN (0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99) DO  (..\sim.exe %lambda% %v0% %r0% %tw% %simDur% %dT% a e 1 %%a %stageOrTime% %log_stages% 0)

