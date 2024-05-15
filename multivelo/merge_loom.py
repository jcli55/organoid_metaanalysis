import loompy

loom_files = [
"/storage/chentemp/u250758/organoid_metaanalysis/Chen/Multi_H9_D35/Multi_H9_D35.loom",
"/storage/chentemp/u250758/organoid_metaanalysis/Chen/Multi_H9_D47/Multi_H9_D47.loom",
"/storage/chentemp/u250758/organoid_metaanalysis/Chen/Multi_organoid_NRL_D075/Multi_organoid_NRL_D075.loom",
"/storage/chentemp/u250758/organoid_metaanalysis/Chen/Multi_H9_D100/Multi_H9_D100.loom",
"/storage/chentemp/u250758/organoid_metaanalysis/Chen/Multi_NRL_GFP_D123/Multi_NRL_GFP_D123.loom",
"/storage/chentemp/u250758/organoid_metaanalysis/Chen/Multi_organoid_D133/Multi_organoid_D133.loom",
"/storage/chentemp/u250758/organoid_metaanalysis/Chen/Multi_Organoid_D183/Multi_organoid_D183.loom",
"/storage/chentemp/u250758/organoid_metaanalysis/Chen/Multi_organoid_NRL_D206/Multi_organoid_NRL_D206.loom",
"/storage/chentemp/u250758/organoid_metaanalysis/Chen/Multi_organoid_NRL_D243/Multi_organoid_NRL_D243.loom"
]
loompy.combine(
    loom_files,
    "/storage/chentemp/u250758/organoid_metaanalysis/Chen/merged.loom",
)
