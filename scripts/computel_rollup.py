import pandas
import numpy
import os


def test():
    print("test")


def open_output(path, sep):
    df = pandas.read_csv(path, sep=sep)
    return(df)


if __name__=="__main__":
    sample_directories = os.listdir('computel')
    computel_columns = ['Sample',
                        'coverage.file',
                        'reads',
                        'read.length',
                        'pattern.length',
                        'base.cov',
                        'num.haploid.chr',
                        'tel.length',
                        'genome.length',
                        'min.seed']

    outputlist = []
    outputlist.append(computel_columns)

    for sample in sample_directories:
        row = [sample.replace('.txt','')]
        path_to_file = 'computel/' + sample + '/tel.length.xls'
        with(open(path_to_file, 'r')) as output_file:
                for line in output_file.readlines():
                    row.append(line.split('\t')[1].replace('\n',''))
        outputlist.append(row)

    print(outputlist)
    with open('computel_summary/computel_rollup.txt', 'w') as rollup_file:
        for row in outputlist:
            rollup_file.write('\t'.join(row) + '\n')
