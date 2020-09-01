import numpy
import pandas
import os
import matplotlib


def test():
    print("test")


def open_output(path, sep):
    df = pandas.read_csv(path, sep=sep)
    return(df)


if __name__=="__main__":
    matplotlib.use('Agg')
    telseq = open_output('telseq_merge/telseq_summary.tsv','\t')
    print(telseq.columns.values)
    print(os.getcwd())
    telomerecat = open_output('telomerecat_merge/telomerecat_summary.csv', ',')
    telomerecat['Sample'] = telomerecat['Sample'].str.replace('.bam','')
    telomerecat = telomerecat.convert_objects(convert_numeric = True)

    computel = open_output('computel_summary/computel_rollup.txt', '\t')

    rollup = pandas.merge(telseq
                         ,telomerecat
                         , on = 'Sample'
                         )

    rollup = pandas.merge(rollup
                         , computel
                         , on = 'Sample'
                         )

    rollup.to_csv('rollup.tsv', sep='\t')

    ax1 = rollup.plot.scatter(x='LENGTH_ESTIMATE',
                             y='Length',
                             c='DarkBlue')

    ax1.get_figure().savefig('telomere_summary/telseq_telomerecat_scatter.pdf')
