To update Figures container:

mcc -m generateFigures.m \
    -a redblue.m \
    -a expfit_delta.m \
    -a eigenfunction_validation.m

To create matlab runtime in apptainer:

apptainer build matlab-runtime-r2025b.sif \
    docker://containers.mathworks.com/matlab-runtime:r2025b
    
Copy to cluster: collectFigures.sh; .sif; -m executable; -m shell

move data to separate dirs:
mkdir ../osgX
mv _* ../osgX/

copy to multiple dirs:
for dir in osg*; do
    cp base/collectFigures.sh base/generateFigures base/run_generateFigures.sh base/collectFigures.sub base/matlab-runtime-r2025b.sif "$dir/"
done

Check dir size (<5GB) before running:
du -h -d 1 | sort -hr

Run the following from /ospool/ap40/data/jovan.zigic:
for dir in osg*/; do
    (cd "$dir" && condor_submit collectFigures.sub)
done

Check status:
condor_q

To cancel:
condor_rm [JOB ID}

Rename after all jobs are done to mark completion:
for dir in osg*; 
    do mv "$dir" "p$dir" 
done

