import sqlite3
import pandas as pd
import py2neo
from py2neo import Graph, Node, Relationship, RelationshipMatcher, NodeMatcher
conn = sqlite3.connect('../resources/STpair.db')
cursor = conn.cursor()

# fetch attributes of datasets
cursor.execute("PRAGMA table_info(datasets)")
attributes = [row[1] for row in cursor.fetchall()]

# fetch all data
cursor.execute("SELECT * FROM datasets")
rows = cursor.fetchall()
cursor.close()
conn.close()

df = pd.DataFrame(rows, columns=attributes)
graph = Graph("bolt://localhost:7687",auth=("neo4j", "biorzh123456"))
df_sc = df[[True if df['dataset_id'][i].startswith("SC") else False for i in df.index]]
df_st = df[[True if df['dataset_id'][i].startswith("ST") else False for i in df.index]]

# create nodes
nodes_st = [Node("ST", id=df_st['dataset_id'][i], title=df_st['title'][i], tissues=df_st['tissues'][i], species=df_st['species'][i], technologies=df_st['technologies'][i]) for i in df_st.index]
nodes_sc = [Node("SC", id=df_sc['dataset_id'][i], title=df_sc['title'][i], tissues=df_sc['tissues'][i], species=df_sc['species'][i], technologies=df_sc['technologies'][i]) for i in df_sc.index]
for node in nodes_st:
  graph.create(node)
for node in nodes_sc:
  graph.create(node)

# create relations
df_st.index = df_st['dataset_id']
df_sc.index = df_sc['dataset_id']
nodes_sc = pd.Series(nodes_sc, index=df_sc['dataset_id'])
nodes_st = pd.Series(nodes_st, index=df_st['dataset_id'])
has_scpaired = df_st['has_paired'].isin(df_sc.index)
has_stpaired = df_st['has_paired'].isin(df_st.index)
relations = []
for id in df_st.index:
  if has_scpaired[id]:
    relations.append(Relationship(nodes_st[id], "PAIR", nodes_sc[df_st['has_paired'][id]]))
  if has_stpaired[id]:
    relations.append(Relationship(nodes_st[id], "PAIR", nodes_st[df_st['has_paired'][id]]))

for relation in relations:
  graph.create(relation)