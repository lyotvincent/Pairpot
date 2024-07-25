from lxml import etree 
import requests
from datetime import datetime
import re

def queryGSE(GSECode, tech='10x Visium'):
  url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GSECode}"  # 网页地址
  response = requests.get(url)  # 发起GET请求
  html = response.text
  lxml_tree = etree.HTML(html)
  record = {}

  # find easy records
  items = ['Title', "Summary", 'Overall design']
  sql_items = ['title', 'summary', 'overall_design']
  for i, item in enumerate(items):
    easy_text = lxml_tree.xpath(f'//tr[td[text()="{item}"]]/td[2]/text()')[0].strip() 
    record[sql_items[i]] = easy_text

  # find organisms
  organisms_text = lxml_tree.xpath('//tr[td[text()="Organism"]]/td[2]/a/text()')
  if len(organisms_text) > 0:
    record['species'] = organisms_text[0]

  # find contributors
  con_text = [a.strip() for a in lxml_tree.xpath('//tr[td[text()="Contributor(s)"]]/td[2]/a/text()')]
  con_text = ','.join(con_text)
  record['contributors'] = con_text

  # find citations
  pmid_text = lxml_tree.xpath('//tr[td[text()="Citation(s)"]]/td[2]/span/a/text()')
  if len(pmid_text) > 0:
    record['pmid'] = pmid_text[0]

  # accessions
  record['accessions'] = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GSECode}"

  # find platforms
  for td in lxml_tree.xpath('//td'):  
      if td.text and 'Platforms' in td.text:  
          platforms_td = td  
          break
  platforms_text = platforms_td.xpath('./following-sibling::td/table/tr/td[2]/text()')
  if len(platforms_text) > 0:
    platforms_text = platforms_text[0].split(' (')[0]
    record['platforms'] = platforms_text

  # find samples
  for td in lxml_tree.xpath('//td'):  
      if td.text and 'Samples' in td.text:  
          samples_td = td  
          break
  samples_text = samples_td.text
  
  samples_text = re.findall("\((\d+)\)", samples_text)
  if len(samples_text) > 0:
    record['n_samples'] = samples_text[0]

  # find dates
  _submission_date = lxml_tree.xpath('//tr[td[text()="Submission date"]]/td[2]/text()')[0]
  parsed_date = datetime.strptime(_submission_date, '%b %d, %Y')
  submission_date = parsed_date.strftime('%Y-%m-%d 00:00:00') 
  _last_update_date = lxml_tree.xpath('//tr[td[text()="Last update date"]]/td[2]/text()')[0]
  parsed_date = datetime.strptime(_last_update_date , '%b %d, %Y')
  last_update_date = parsed_date.strftime('%Y-%m-%d 00:00:00') 
  record['submission_date'] = submission_date
  record['last_modified'] = last_update_date

  # find contacts
  contact_tr = lxml_tree.xpath('//tr[td[text()="E-mail(s)"]]/td[2]/a/text()')
  if len(contact_tr) > 0:
    record['contacts'] = contact_tr[0]

  # default attr
  record['cells'] = 0
  record['spots'] = 0
  record['genes'] = 0
  record['sex'] = "Unknown"
  record['technologies'] = tech
  record

  sql_cols = ["id", "dataset_id", "title", "species", "tissues", "organ_parts", "cell_types", "cells", "spots", "genes", "development_stages", "sex", "technologies", "n_samples", "n_sections", "disease", "summary", "overall_design", "submission_date", "last_modified", "contributors", "contacts", "citation", "accessions", "platforms", "pmid", "has_paired"]
  values = ['']*len(sql_cols)
  for i,col in enumerate(sql_cols):
      if col in record.keys():
          values[i] = record[col]
  values_str = ','.join([f"\'{str(col)}\'" for col in values])
  sql_str = f"INSERT INTO datasets ({','.join([col for col in sql_cols])}) VALUES ({values_str})"
  return record, sql_str