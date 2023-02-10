if [ "$SRC_SURFEX" == "" ] 
then
	echo "Problem: SRC_SURFEX needs to be defined. Run the appropriate profile-surfex... file."
elif [ "$OPTLEVEL" == "" ]
then
	echo "Problem: OPTLEVEL needs to be defined. Run the appropriate profile-surfex... file."
elif [ "$VER_MPI" == "" ]
then
 	echo "Problem: VER_MPI needs to be defined. Run the appropriate profile-surfex... file."
else
#removes old links
	rm -f *.exe
	rm -f *.bin
	rm -f \**

#links executables
  	ln -s $SRC_SURFEX/exe/PGD-LXgfortran-SFX-V8-1-1-$VER_MPI-$VER_OMP-$OPTLEVEL-X0 pgd.exe
  	ln -s $SRC_SURFEX/exe/PREP-LXgfortran-SFX-V8-1-1-$VER_MPI-$VER_OMP-$OPTLEVEL-X0 prep.exe
	ln -s $SRC_SURFEX/exe/OFFLINE-LXgfortran-SFX-V8-1-1-$VER_MPI-$VER_OMP-$OPTLEVEL-X0 offline.exe
	ln -s $SRC_SURFEX/exe/SODA-LXgfortran-SFX-V8-1-1-$VER_MPI-$VER_OMP-$OPTLEVEL-X0 soda.exe
	#ln -s $SRC_SURFEX/exe/SXPOST-LXgfortran-SFX-V8-1-1-$VER_MPI-$VER_OMP-$OPTLEVEL-X0 sxpost.exe

#links ecoclimap .bin files
	ln -s $SRC_SURFEX/MY_RUN/ECOCLIMAP/*.bin .

#links physiography files
        ln -s $HOME/SFX/data/* .
fi
