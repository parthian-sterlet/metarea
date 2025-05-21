#!/bin/sh

g++ -o pauc_forback_2motifs.exe pauc_forback_2motifs.cpp
g++ -o pauc_forback_2motif0.exe pauc_forback_2motif0.cpp
g++ -o pauc_forback_2motif0.exe pauc_forback_2motif0s.cpp
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
chmod a+x pauc_forback_2motif0.exe
chmod a+x pauc_forback_2motif0s.exe
chmod a+x sitega_thr_dist_mat.exe
chmod a+x pfm_to_pwm_mat.exe
chmod a+x pwm_iz_pwm_thr_dist0.exe

cd ../run
chmod a+x com_line_pwm_sga
chmod a+x com_line_pwm_pwm
chmod a+x com_line_anc_lib
chmod a+x com_line_lib_lib
chmod a+x com_line_de_novo_pwm
chmod a+x com_line_de_novo_sga
chmod a+x den_pwm.pl
chmod a+x den_sga.pl

cd ../partners
cat h12core_hg38.binary.tar.gz.part* > h12core_hg38.binary.tar.gz
cat h12core_mm10.binary.tar.gz.part* > h12core_mm10.binary.tar.gz
tar -xvzf h12core_hg38.binary.tar.gz
tar -xvzf h12core_mm10.binary.tar.gz

cd ../genomes/at
tar -xvzf ups1500_at10.seq.tar.gz
cd ../dm
tar -xvzf ups1500_dm6.seq.tar.gz
cd ../mm
tar -xvzf ups2kb_mm10.seq.tar.gz
cd ../hs
tar -xvzf ups2kb_hg38.seq.tar.gz
cd ..
cd ..
