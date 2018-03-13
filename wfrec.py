import intervaltree as it
import msprime
import numpy as np
import sys


def wfrec(nsam, rho, nsites, theta):
    samples = []
    for i in range(nsam):
        samples.append(it.IntervalTree([it.Interval(0, nsites)]))

    links = np.array(
        [i[1] - i[0] - 1 for j in samples for i in j], dtype=np.int)
    nlinks = links.sum()

    n = nsam
    rbp = rho / float(nsites - 1)
    t = 0.0

    nodes = msprime.NodeTable()
    edges = msprime.EdgeTable()

    nodes.set_columns(time=np.zeros(
        nsam), flags=np.ones(nsam, dtype=np.uint32))

    sample_indexes = [i for i in range(len(samples))]
    next_index = len(sample_indexes)

    while(n > 1):
        rcoal = float(n * (n - 1))
        rrec = rbp * float(nlinks)

        iscoal = bool(np.random.random_sample(1)[0] < rcoal / (rcoal + rrec))
        t += np.random.exponential(4. / (rcoal + rrec), 1)[0]
        assert len(samples) == len(links), "sample/link error"
        if iscoal is True:
            chroms = np.sort(np.random.choice(n, 2, replace=False))
            c1 = chroms[0]
            c2 = chroms[1]

            nodes.add_row(time=t, flags=msprime.NODE_IS_SAMPLE)
            for i in samples[c1]:
                edges.add_row(left=i[0], right=i[1],
                              parent=next_index, child=sample_indexes[c1])
                edges.add_row(left=i[0], right=i[1],
                              parent=next_index, child=sample_indexes[c2])
            newchrom = it.IntervalTree()
            # Merge intervals of the two chromosomes
            # and remove overlaps
            for i in samples[c1]:
                newchrom.append(i)
            for i in samples[c2]:
                newchrom.append(i)
            newchrom.merge_overlaps()
            samples.pop(c2)
            samples.pop(c1)
            samples.append(newchrom)
            sample_indexes.pop(c2)
            sample_indexes.pop(c1)
            sample_indexes.append(next_index)
            next_index += 1
            n -= 1
        else:
            # Pick a chrom proportional to
            # its total size:
            print("check", len(samples), len(sample_indexes), len(links))
            chrom = np.random.choice(
                len(sample_indexes), 1, p=links / links.sum())[0]
            print(chrom)
            print(len(samples[chrom]))
            mnpos = min([i for j in samples[chrom]
                         for i in j if i is not None])
            mxpos = max([i for j in samples[chrom]
                         for i in j if i is not None])
            pos = np.random.randint(mnpos, mxpos)
            print("to chop",samples[chrom])
            samples[chrom].chop(pos, pos)
            print("chopped",samples[chrom])
            tc = it.IntervalTree([i for i in samples[chrom] if i[0] >= pos])
            samples[chrom].remove_overlap(pos, nsites)
            if len(samples[chrom]) == 0:
                print(pos,tc)
                sys.exit(0)
            samples.append(tc)
            sample_indexes.append(next_index)
            print(samples)
            next_index += 1
            n += 1

        assert len(samples) == len(sample_indexes), "sample/sample_index error"
        links = np.array(
            [i[1] - i[0] - 1 for j in samples for i in j], dtype=np.int)
        nlinks = links.sum()
        # print(samples)
        # print(len(samples),len(links))
        assert len(samples) == len(links), "sample/link error 2"

    msprime.sort_tables(nodes=nodes, edges=edges)
    return msprime.load_tables(nodes=nodes, edges=edges)


def test():
    S = []
    np.random.seed(42)
    msp_rng = msprime.RandomGenerator(84)
    for i in range(1000):
        ts = wfrec(100, 0, 1000, 100)
        sites = msprime.SiteTable()

        mutations = msprime.MutationTable()
        mutgen = msprime.MutationGenerator(
            msp_rng, 100. / 4000.)
        mutgen.generate(ts.tables.nodes, ts.tables.edges, sites, mutations)
        ts = ts.load_tables(nodes=ts.tables.nodes,
                            edges=ts.tables.edges,
                            sites=sites,
                            mutations=mutations)
        S.append(ts.num_mutations)
    S2 = []
    for i in msprime.simulate(100, recombination_rate=0, mutation_rate=25, num_replicates=1000):
        # S2.append(i.tables.nodes[next(i.trees()).root].time)
        S2.append(i.num_mutations)
    return S, S2
