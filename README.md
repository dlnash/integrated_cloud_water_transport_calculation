## ICT, IWT, and ICT/IWT Plots

This repository uses GFS data to compute and plot integrated condensate transport (ICT), integrated w? transport (IWT), and the Ratio of ICT to integrated water vapor transport (IVT) for the Interior West. 

### To run:

---

To run the code with a singularity container:

```bash
## runs plots for GFS
singularity exec --bind /data:/data,/home:/home,/work:/work,/common:/common -e /data/projects/operations/ICT_IWT/envs/ICT_IWT.sif /opt/conda/bin/python /data/projects/operations/ICT_IWT/run_tool.py "GFS" "YYYYMMDDHH"
```