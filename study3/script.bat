@echo off

set lambda=1
set v0=1
set r0=1
set simDur=2e6
set dT=0.01
set stageOrTime=1
set log_stages=0

REM tw=pi/16
START FOR %%a IN (0.01 0.5 1 1.5 2 2.5 3 3.5 4) DO  (..\sim.exe %lambda% %v0% %r0% 0.1963 %simDur% %dT% a d %%a 0 %stageOrTime% %log_stages% 0)

REM tw=2*pi/16
START FOR %%a IN (0.01 0.5 1 1.5 2 2.5 3 3.5 4) DO  (..\sim.exe %lambda% %v0% %r0% 0.3927 %simDur% %dT% a d %%a 0 %stageOrTime% %log_stages% 0)

REM tw=3*pi/16
START FOR %%a IN (0.01 0.5 1 1.5 2 2.5 3 3.5 4) DO  (..\sim.exe %lambda% %v0% %r0% 0.5890 %simDur% %dT% a d %%a 0 %stageOrTime% %log_stages% 0)

REM tw=4*pi/16
START FOR %%a IN (0.01 0.5 1 1.5 2 2.5 3 3.5 4) DO  (..\sim.exe %lambda% %v0% %r0% 0.7854 %simDur% %dT% a d %%a 0 %stageOrTime% %log_stages% 0)

REM tw=5*pi/16
START FOR %%a IN (0.01 0.5 1 1.5 2 2.5 3 3.5 4) DO  (..\sim.exe %lambda% %v0% %r0% 0.9817 %simDur% %dT% a d %%a 0 %stageOrTime% %log_stages% 0)

REM tw=6*pi/16
START FOR %%a IN (0.01 0.5 1 1.5 2 2.5 3 3.5 4) DO  (..\sim.exe %lambda% %v0% %r0% 1.1781 %simDur% %dT% a d %%a 0 %stageOrTime% %log_stages% 0)

REM tw=7*pi/16
START FOR %%a IN (0.01 0.5 1 1.5 2 2.5 3 3.5 4) DO  (..\sim.exe %lambda% %v0% %r0% 1.3744 %simDur% %dT% a d %%a 0 %stageOrTime% %log_stages% 0)

REM tw=8*pi/16
START FOR %%a IN (0.01 0.5 1 1.5 2 2.5 3 3.5 4) DO  (..\sim.exe %lambda% %v0% %r0% 1.5708 %simDur% %dT% a d %%a 0 %stageOrTime% %log_stages% 0)