SHELL := /bin/bash

# all class 1, 2 runs
CLASS12_ALL_RUNS = class12 class12_tpcf_nbins class12_crossf_nbins
# all class 1, 2, 3 runs
CLASS123_ALL_RUNS = class123 class123_tpcf_nbins class123_crossf_nbins

class12_all:  $(CLASS12_ALL_RUNS)
class123_all: $(CLASS123_ALL_RUNS)

run01: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'run01'" sc_class="'class12'" do_nn="True" > run01.log 2>&1
	cp run01.log run01_`date +%Y-%m-%dT%H:%M`.log

class12: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'class12'" sc_class="'class12'" do_catalog_setup="True" do_nn="True" do_assoc="True" do_tpcf="False" > class12.log 2>&1
	cp class12.log class12_`date +%Y-%m-%dT%H:%M`.log

class12_nnfix: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'class12'" sc_class="'class12'" nbins_crosscf="[5,10,15,20,25]" do_catalog_setup="False" do_nn="True" do_assoc="False" do_tpcf="False" do_cf="False"> class12.log 2>&1
	cp class12.log class12_`date +%Y-%m-%dT%H:%M`.log

class12_assocfix: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'class12'" sc_class="'class12'" nbins_crosscf="[5,10,15,20,25]" do_catalog_setup="False" do_nn="False" do_assoc="True" do_tpcf="False" do_cf="False"> class12.log 2>&1
	cp class12.log class12_`date +%Y-%m-%dT%H:%M`.log

# run after running 'class12' to play with number of bins for the two point correlation functions
class12_tpcf_nbins: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'class12'" sc_class="'class12'" nbins_tpcf="[5,10,15,20,25]" do_catalog_setup="False" do_nn="False" do_assoc="False" do_tpcf="True" > class12.log 2>&1
	cp class12.log class12_`date +%Y-%m-%dT%H:%M`.log

# run after running 'class12' to play with number of bins for the two point correlation functionsG
class12_crossf_nbins: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'class12'" sc_class="'class12'" nbins_crosscf="[5,10,15,20,25]" do_catalog_setup="False" do_nn="False" do_assoc="False" do_tpcf="False" do_cf="True"> class12.log 2>&1
	cp class12.log class12_`date +%Y-%m-%dT%H:%M`.log

# with class 1, 2, 3 clusters
class123: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'class123'" sc_class="'class123'" do_catalog_setup="True" do_nn="True" do_assoc="True" do_tpcf="False" > class123.log 2>&1
	cp class123.log class123_`date +%Y-%m-%dT%H:%M`.log

class123_nnfix: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'class123'" sc_class="'class123'" nbins_crosscf="[5,10,15,20,25]" do_catalog_setup="False" do_nn="True" do_assoc="False" do_tpcf="False" do_cf="False"> class123.log 2>&1
	cp class123.log class123_`date +%Y-%m-%dT%H:%M`.log

class123_assocfix: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'class123'" sc_class="'class123'" nbins_crosscf="[5,10,15,20,25]" do_catalog_setup="False" do_nn="False" do_assoc="True" do_tpcf="False" do_cf="False"> class123.log 2>&1
	cp class123.log class123_`date +%Y-%m-%dT%H:%M`.log

class123_tpcf_nbins: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'class123'" sc_class="'class123'" nbins_tpcf="[5,10,15,20,25]" do_catalog_setup="False" do_nn="False" do_assoc="False" do_tpcf="True" > class123.log 2>&1
	cp class123.log class123_`date +%Y-%m-%dT%H:%M`.log

# run after running 'class12' to play with number of bins for the two point correlation functionsG
class123_crossf_nbins: 
	time python full_run.py galaxy_list_file="'galaxy_go.list'" run_name="'class123'" sc_class="'class123'" nbins_crosscf="[5,10,15,20,25]" do_catalog_setup="False" do_nn="False" do_assoc="False" do_tpcf="False" do_cf="True"> class123.log 2>&1
	cp class123.log class123_`date +%Y-%m-%dT%H:%M`.log
