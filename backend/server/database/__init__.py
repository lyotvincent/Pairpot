import sqlite3
dataset_columns = ['id', 'dataset_id', 'title', 'species', 'tissues', 'organ_parts',
                   'cell_types', 'cells', 'spots', 'genes', 'development_stages', 'sex',
                   'technologies', 'n_samples', 'n_sections', 'disease', 'summary',
                   'overall_design', 'submission_date', 'last_modified', 'contributors',
                   'contacts', 'citation', 'accessions', 'platforms', 'pmid', 'has_paired']



# 连接到SQLite数据库
conn = sqlite3.connect('resources/STpair.db')

# 创建一个游标对象
cursor = conn.cursor()

# 执行PRAGMA命令来获取表的属性名
cursor.execute("PRAGMA table_info(datasets)")

# 获取结果并输出属性名
attributes = [row[1] for row in cursor.fetchall()]
print(attributes)

conn.close()