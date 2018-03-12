import intervaltree as it
import msprime
import numpy as np


def wfrec(nsam, rho, nsites, theta):
    samples = []
    for i in range(nsam):
        samples.append(it.IntervalTree([it.Interval(0, nsites)]))

    nlinks = (nsites - 1) * nsam
    # nlinks2 = sum([i[1]-i[0]-1 for j in samples for i in j])

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

        tcoal = np.random.exponential(4. / rcoal, 1)[0]

        t += tcoal

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

    msprime.sort_tables(nodes=nodes, edges=edges)
    return msprime.load_tables(nodes=nodes, edges=edges)


def test():
    S = []
    np.random.seed(42)
    msp_rng = msprime.RandomGenerator(84)
    for i in range(1000):
        ts = wfrec(100, 100, 1000, 100)
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
    for i in msprime.simulate(100, mutation_rate=25, num_replicates=1000):
        # S2.append(i.tables.nodes[next(i.trees()).root].time)
        S2.append(i.num_mutations)
    return S, S2
