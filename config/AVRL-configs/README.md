# NextFlow Configuration Files to Make Switching Between Computers Easy

When bundled together with your nextflow workflow of choice, these config files can configure scripts for execution on CHTC, dhogal(2), or an iMac desktop of either available architecture. Simply clone or copy the config file associated with the computer you're using and add it into your nextflow run command, like so:

```
nextflow run workflow.nf -c /absolute/or/relative/path/to/chtc.config
```