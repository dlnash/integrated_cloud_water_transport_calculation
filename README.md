# Repository Summary

This repository uses GFS data to compute and visualize water transport metrics for the Interior West, including:

- **Integrated Water Vapor Transport (IVT):** Transport of atmospheric water vapor.  
- **Integrated Water Transport (IWT):** IVT calculated using total mixing ratios, including condensate.  
- **Integrated Condensate Transport (ICT):** The difference between IWT and IVT.  
- **ICT/IVT Ratio (%):** Percentage of condensate transport relative to vapor transport.


## To run:

---

To run the code with a singularity container:

```bash
## runs plots for GFS
apptainer exec -e --bind /data:/data /data/projects/operations/ICT_IWT/env/ICT_IWT.sif python /data/projects/operations/ICT_IWT/run_tool.py "GFS" "YYYYMMDDHH"
```