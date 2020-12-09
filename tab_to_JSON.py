#!/usr/bin/env python
import sys
import json
from collections import defaultdict

#Input results file
Kinase_infile = open("output/Annotated_Merge_top5.csv", 'r')

#JSON outfile *In configurations*
JSON_outfile = open(sys.argv[1], 'w')

data = []
edges = []
kinases = dict()


def addNode(node_type, node):
    classes = [node_type]
    if node in kinases:
        classes.append("kinase")
        if kinases[node]: classes.append("dark")

    json_node = {
        "data": {
            "id": node
        },
        "classes": classes,
    }
    data.append(json_node)


def addEdge(source, target, edge_type, i, numBaits):
    json_edge = {
        "data": {
            "id": "e" + str(i),
            "source": source,
            "target": target,
            "numBaits": numBaits
        },
        "classes": [edge_type]
    }
    data.append(json_edge)


for line in Kinase_infile:
    pl = line.strip().split(",")
    prey_prot_for_kinase = pl[4]
    if not "NA" in pl[50]:
        is_dark = pl[50] == "Dark"
        kinases[prey_prot_for_kinase] = is_dark

baits = set()
prey = defaultdict(dict)

Kinase_infile.seek(0)

header = Kinase_infile.readline().split(',')
for line in Kinase_infile:
    pl = line.split(",")
    bait_gene = pl[58]
    AP_type = "MiniTurbo"
    prey_prot = pl[4]
    prey_gene = pl[14]

    if bait_gene != prey_gene and prey_gene != '':
        baits.add(bait_gene)
        if not AP_type in prey[prey_gene]:
            prey[prey_gene][AP_type] = 1
        else:
            prey[prey_gene][AP_type] += 1

        edges.append([bait_gene, prey_gene, AP_type])

prey_set = set(prey.keys()) - baits

for node in baits: addNode("bait", node)
for node in prey: addNode("prey", node)

for i in range(0, len(edges)): addEdge(edges[i][0], edges[i][1], edges[i][2], i, prey[edges[i][1]][edges[i][2]])

# Write JSON
JSON_data = json.dumps(data)
JSON_outfile.write(JSON_data)
