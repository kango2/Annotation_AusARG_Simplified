# Annotation_AusARG_Simplified
This repository is basically a cleaner version of [this repo](https://github.com/kango2/Annotation_AusARG) and contains a simple script to automatically launch the pipeline. Everything important can be set in this file:
```
launch.sh
```

# Usage
To get started, clone this repository to somewhere you have access.
```
cd /path/to/somewhere/you/like
git clone https://github.com/kango2/Annotation_AusARG_Simplified.git
```
Then, open `/path/to/somewhere/you/like/Annotation_AusARG_Simplified/launch.sh` in any text editor and configure your settings (any lines starting with `export` in the first 43 lines), Once you are done with this and saved the file, launch it.
```
cd /path/to/somewhere/you/like/Annotation_AusARG_Simplified
./launch.sh
```
That's it, the script should have automatically launched 2 or 3 PBS jobs for you (depends on whether you are running P2 or P1). You can check status of these job in the queue with
```
qstat
qstat -u ${USER}
```
