#!/usr/bin/env python3
import Groups
import runAlignment

def buildGroupAlignments(D, namD, rev_namD, F, pD):
    groups = Groups.groupAlignmentDict(D, namD)
    aliD = dict()
    for i, group in enumerate(groups):
        subF = dict((nam, F[rev_namD[nam]]) for nam in group)
        ali = runAlignment.runAlignmentPW(subF, pD, quick=False,
                                          keep_failed=True)
        aliD[i] = ali
    return (aliD)