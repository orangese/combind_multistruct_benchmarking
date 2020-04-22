import click
import pandas as pd


@click.command()
@click.option('--n', default=10)
@click.option('--cut', default=100)
@click.argument('input_file')
@click.argument('output_file')
def main(input_file, output_file, n, cut):
	df = pd.read_csv(input_file)
	df = df.loc[df['AFFINITY'] < cut]
	df = df.sample(n)

	df.to_csv(output_file, index=False)

main()
