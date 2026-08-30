To update Figures container:

mcc -m generateFigures.m \
    -a redblue.m \
    -a expfit_delta.m \
    -a eigenfunction_validation.m

To create matlab runtime in apptainer:

apptainer build matlab-runtime-r2025b.sif \
    docker://containers.mathworks.com/matlab-runtime:r2025b
    
Copy to cluster: .sif; -m executable; -m shell

