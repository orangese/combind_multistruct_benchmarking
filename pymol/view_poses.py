from pymol import cmd
import sys
print sys.argv

def load_results(protein):
    cmd.load("glide_{}_pv.mae".format(protein))
    cmd.load("combind_{}_pv.mae".format(protein))
    cmd.load("best_{}_pv.mae".format(protein))
    cmd.load("crystal_{}_pv.mae".format(protein))
    
    cmd.util.cbay("glide_{}_pv".format(protein))
    cmd.util.cbac("combind_{}_pv".format(protein))
    cmd.util.cbam("best_{}_pv".format(protein))
    cmd.util.cbag("crystal_{}_pv".format(protein))
    cmd.util.cbab("not het")

def show_prot():
    cmd.enable("*_prot")

def hide_prot():
    cmd.disable("*_prot")

def show_lig(lig):
    cmd.disable('*.*')
    cmd.enable("*{}_lig".format(lig))
    cmd.zoom("*{}_lig".format(lig))

cmd.extend('load_results', load_results)
cmd.extend('show_lig', show_lig)
cmd.extend('show_prot', show_prot)
cmd.extend('hide_prot', hide_prot)
