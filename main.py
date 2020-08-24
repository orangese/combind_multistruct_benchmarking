#!/bin/env run
import os
import click
import utils
import config

@click.group()
@click.option('--data', default='/oak/stanford/groups/rondror/users/jpaggi/combind')
@click.option('--ligands', default='{ROOT}/structures/pdb.csv')
@click.pass_context
def main(ctx, data, ligands):
    paths = {'CODE': os.path.dirname(os.path.realpath(__file__)),
             'DATA': data,
             'PDB': ligands}
    paths.update(config.PATHS)
    paths = utils.resolve(paths)
    ctx.obj = paths

@main.command()
@click.option('--stats-version', default='rd1')
@click.option('--struct', default=None)
@click.option('--processes', default=1)
@click.option('--more-ligands', default=None)
@click.argument('task')
@click.argument('proteins', nargs=-1)
@click.pass_obj
def prepare(paths, task, stats_version, proteins, struct, processes, more_ligands):
    """
    Prep structures and ligands; docking; and featurization.
    """
    import dock.prepare_all
    params = config.STATS[stats_version]
    proteins = list(proteins)
    if not proteins:
        proteins = utils.get_proteins(paths, [])
    dock.prepare_all.main(params, paths, task, list(proteins), struct, processes, more_ligands)

@main.command()
@click.option('--stats-version', default='rd1')
@click.option('--stats-root', default='{}/stats_data/rd1'.format(os.environ['COMBINDHOME']))
@click.option('--struct', default=None)
@click.option('--pose-fname', default='poses.sc')
@click.option('--alpha', default=1.0)
@click.option('--gc50', default=10.0)
@click.option('--num-poses', default=100)
@click.option('--features', default='mcss,hbond,sb,contact')
@click.option('--max-iterations', default=int(1e6))
@click.option('--restart', default=500)
@click.option('--xtal', multiple=True)
@click.argument('protein')
@click.argument('queries', nargs=-1)
@click.pass_obj
def score(paths, stats_root, protein, struct, queries, pose_fname, xtal,
          stats_version, alpha, gc50, num_poses, features, max_iterations, restart):
    """
    Run ComBind!
    """
    import score.scores
    queries = list(queries)
    xtal = list(xtal)
    params = config.STATS[stats_version]
    features = {feature: config.FEATURE_DEFS[feature]
                for feature in features.split(',')}

    score.scores.score(paths, params, features, stats_root, protein, queries,
                       xtal=xtal, struct=struct, pose_fname=pose_fname,
                       alpha=alpha, gc50=gc50, num_poses=num_poses,
                       max_iterations=max_iterations, restart=restart)

@main.command()
@click.option('--merged-root')
@click.option('--stats-version', default='default')
@click.argument('stats-root')
@click.argument('proteins', nargs=-1)
@click.pass_obj
def statistics(paths, stats_version, stats_root, proteins, merged_root):
    """
    Fit statistics to a set of protein–ligand complexes.
    """
    import score.statistics
    params = config.STATS[stats_version]
    proteins = list(proteins)
    if not proteins:
        proteins = utils.get_proteins(paths, [])
    score.statistics.compute(params, paths, config.FEATURE_DEFS, stats_root,
                             proteins, merged_root)

@main.command()
@click.argument('merged')
def plot_statistics(merged):
    import score.statistics
    score.statistics.plot(merged)

@main.command()
@click.argument('ifp_version')
@click.argument('input_file')
@click.argument('output_file')
@click.argument('poses', type=int)
def ifp(ifp_version, input_file, output_file, poses):
    """
    Compute interaction fingerprints. (For internal use.)
    """
    import ifp.ifp
    ifp.ifp.ifp(config.IFP[ifp_version], input_file, output_file, poses)

@main.command()
@click.argument('mode')
@click.argument('ligand1')
@click.argument('ligand2')
@click.argument('ligand1_path')
@click.argument('ligand2_path')
@click.argument('init_file')
@click.argument('mcss_types_file')
@click.argument('rmsd_file', default='')
@click.argument('max_poses', type=int, default=0)
@click.argument('mcss_string_rep', default='')
def mcss(mode, ligand1, ligand2, ligand1_path, ligand2_path,
         init_file, mcss_types_file, rmsd_file, max_poses,
         mcss_string_rep):
    """
    Compute substructure similiarity. (For internal use.)
    """
    import mcss.mcss
    mcss.mcss.main(mode, ligand1, ligand2, ligand1_path, ligand2_path,
                   init_file, mcss_types_file, rmsd_file, max_poses,
                   mcss_string_rep)

main()
