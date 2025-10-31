import numpy as np

def get_params(array_group):
    ell_values = np.round(np.arange(1.02, 1.03, 0.08), 2)
    params = []
    if array_group == 0:
        K_values = range(0,1)
    elif array_group == 1:
        K_values = []
    else:
        return []  # no other groups

    for argK in K_values:
        if argK == 0:
            argT_values = np.round(np.arange(0.0,2.1,0.1),1)
            argdt,argN,mem=1e-3,24,"32G"
        elif argK == 1:
            argT_values = np.round(np.arange(-0.5,1.6,0.1),1)
            argdt,argN,mem=1e-3,32,"32G"
        elif argK == 2:
            argT_values = np.round(np.arange(-1.0,1.1,0.1),1)
            argdt,argN,mem=1e-3,32,"32G"
        elif argK == 3:
            argT_values = np.round(np.arange(-1.5,0.6,0.1),1)
            argdt,argN,mem=5e-5,48,"32G"
        elif argK == 4:
            argT_values = np.round(np.arange(-2.0,-0.9,0.1),1)
            argdt,argN,mem=1e-6,64,"32G"
        elif argK == 5:
            argT_values = np.round(np.arange(-3.0,-1.9,0.1),1)
            argdt,argN,mem=1e-8,192,"48G"

        for argell in ell_values:
            for argT in argT_values:
                params.append((argK, argell, argT, argdt, argN, mem))
    return params

def total_params(array_group):
    return len(get_params(array_group))

def available_array_groups():
    groups = []
    for g in [0,1]:
        if total_params(g) > 0:
            groups.append(g)
    return groups

if __name__ == "__main__":
    import argparse, os, sys

    parser = argparse.ArgumentParser()
    parser.add_argument("task_id", nargs="?", type=int)
    parser.add_argument("--array-group", type=int, choices=[0,1], default=0)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    params = get_params(args.array_group)

    if args.dry_run:
        print(f"Array group {args.array_group}: total tasks = {len(params)}")
        if len(params) > 0:
            argK,argell,argT,argdt,argN,mem = params[0]  # first task
            jobname = f"run_{argK}_{argell:.2f}_{argT:.1f}_{argdt}_{argN}"
            logfile = f"./output/{jobname}.log"
            matlab_cmd = (
                f"matlab -nodisplay -nodesktop -nosplash -r "
                f"\"main_2DKS({argdt},{argN},{argK},{argK},1,{argell},{argell},0.02,{argT},{argT},1,'optimize','IC',1e-6,0.0); exit\" "
                f"> {logfile} 2>&1 < /dev/null"
            )
            print("Example MATLAB command for first task:")
            print(matlab_cmd)
        sys.exit(0)

    if args.task_id is None:
        print("Error: must provide task_id when not using --dry-run")
        sys.exit(1)

    if args.task_id >= len(params):
        print(f"Error: task_id {args.task_id} out of range")
        sys.exit(1)

    argK,argell,argT,argdt,argN,mem = params[args.task_id]
    jobname = f"run_{argK}_{argell:.2f}_{argT:.1f}_{argdt}_{argN}"
    logfile = f"./output/{jobname}.log"
    os.makedirs("./output", exist_ok=True)

    cmd = (
        f"matlab -nodisplay -nodesktop -nosplash -r "
        f"\"main_2DKS({argdt},{argN},{argK},{argK},1,{argell},{argell},0.02,{argT},{argT},1,'optimize','IC',1e-6,0.0); exit\" "
        f"> {logfile} 2>&1 < /dev/null"
    )
    print(f"Running task {args.task_id} (array group {args.array_group}): {jobname}")
    os.system(cmd)
