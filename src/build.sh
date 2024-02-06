#!/bin/sh

g++ -o pauc_forback_2motifs.exe pauc_forback_2motifs.cpp
g++ -o pauc_forback_anc_lib.exe pauc_forback_anc_lib.cpp
g++ -o pauc_forback_pwm_sga_only.exe pauc_forback_pwm_sga_only.cpp
g++ -o pauc_forback_2motifs_only.exe pauc_forback_2motifs_only.cpp
g++ -o sitega_thr_dist_mat.exe sitega_thr_dist_mat.cpp
g++ -o pfm_to_pwm_mat.exe pfm_to_pwm_mat.cpp
g++ -o pwm_iz_pwm_thr_dist0.exe pwm_iz_pwm_thr_dist0.cpp
chmod a+x pauc_forback_pwm_sga_only.exe
chmod a+x pauc_forback_2motifs_only.exe
chmod a+x pauc_forback_anc_lib.exe
chmod a+x pauc_forback_2motifs.exe
chmod a+x sitega_thr_dist_mat.exe
chmod a+x pfm_to_pwm_mat.exe
chmod a+x pwm_iz_pwm_thr_dist0.exe
chmod a+x run_pwm_sga

cd ..
cd genomes
cd mm10
tar -cvzf ups2kb_mm10.seq.tar.gz ups2kb_mm10.seq
cd ..
cd ..
