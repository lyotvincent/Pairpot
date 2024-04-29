import {
  DatabaseOutlined,
  ExperimentOutlined,
  SwitcherOutlined,
  YuqueOutlined,
  EyeOutlined,
  LinkOutlined,
  FileTextOutlined,
  DownloadOutlined,
} from '@ant-design/icons'
import React, { useRef, useState, useEffect } from 'react'
import {
  Avatar,
  List,
  Space,
  Checkbox,
  Card,
  Tooltip,
  Button,
  Rate,
} from 'antd'
import DatasetDescription from './DatasetDescription'
import TextCollapse from './TextCollapse'
import { useNavigate } from 'react-router-dom'
import axios from 'axios'
//import { response } from 'express'

const IconText = ({ icon, text, attr, onClick, placement = 'bottom' }) => (
  <Tooltip placement={placement} title={attr} arrow={false} align={'center'}>
    <Space onClick={onClick} style={{ cursor: 'pointer' }}>
      {React.createElement(icon)}
      {text}
    </Space>
  </Tooltip>
)

const AddCommas = (num) => {
  return num.toLocaleString().replace(/^0+/, '')
}

const DatasetList = ({ src, col }) => {
  const navigate = useNavigate()
  const [dataSrc, setDataSrc] = useState([])
  const [attr, setAttr] = useState([])
  const [dataCol, setDataCol] = useState({})
  const [descOpen, setDescOpen] = useState(false)
  const [descInfo, setDescInfo] = useState([])
  useEffect(() => {
    setDataSrc(src)
    setDataCol(col)

    var attrStack = {}
    for (var i = 0; i < col.length; i++) {
      var colname = col[i]
      attrStack[colname] = i
    }
    setAttr(attrStack)
  }, [src, col])
  return (
    <List
      itemLayout="vertical"
      size="large"
      pagination={{
        showSizeChanger: true,
        onChange: (page) => {
          console.log(page)
        },
        pageSize: 5,
      }}
      dataSource={dataSrc}
      renderItem={(item) => (
        <List.Item
          key={item[0]}
          actions={[
            <IconText
              icon={YuqueOutlined}
              attr="Species"
              text={item[attr.species] ? item[attr.species] : 'Unknown'}
              key="list-vertical-star-o"
            />,
            <IconText
              icon={SwitcherOutlined}
              attr="Tissue(s)"
              text={item[attr.tissues] ? item[attr.tissues] : 'Unknown'}
              key="list-vertical-like-o"
            />,
            <IconText
              icon={ExperimentOutlined}
              attr="Technology"
              text={item[attr.technologies]}
              key="list-vertical-message"
              onClick={() => {}}
            />,
            <Rate defaultValue={2} />,
          ]}
          extra={
            <div>
              <br />
              <DatasetDescription
                descCol={dataCol}
                descInfo={item}
                text={'Descriptions'}
                placement={'left'}
              />
              <br />
              <IconText
                icon={EyeOutlined}
                text="Visualization"
                key="list-vertical-visualization-o"
                attr="Visualize this Dataset"
                placement="left"
                onClick={()=>{
                  let values = Object.fromEntries(
                    dataCol.map((k, i) => [k, item[i]])
                  )
                  let scid = values['has_paired']
                  let scitem = dataSrc.find((item)=> item[1]===scid)
                  let state = {
                    st: values
                  }
                  if(typeof scitem !== 'undefined'){
                    let scvalues = Object.fromEntries(
                      dataCol.map((k, i) => [k, scitem[i]])
                    )
                    state['sc'] = scvalues
                  }
                  navigate('/browse', { state: state })
                }}
              />
              <br />
              <br />
              <IconText
                icon={LinkOutlined}
                text="Paired Datasets"
                key="list-vertical-Paired-Datasets-o"
                attr="Link paired datasets"
                placement="left"
                onClick={() => {
                  let values = Object.fromEntries(
                    dataCol.map((k, i) => [k, item[i]])
                  )
                  navigate('/submit/link', { state: values })
                  console.log(
                    Object.fromEntries(dataCol.map((k, i) => [k, item[i]]))
                  )
                }}
              />
              <br/>
              <br />
              <IconText
                icon={DownloadOutlined}
                text="Download Datasets"
                key="list-vertical-Download-Datasets-o"
                attr="Download datasets"
                placement="left"
                onClick={() => {
                  axios
                  .get('http://localhost:5522/query', {
                    responseType: 'blob',
                  })
                  .then((response) => {
                    console.log(response)
                    let blob = response.data
                    let url = window.URL.createObjectURL(blob);
                    let a = document.createElement('a');
                    a.href = url;
                    a.download = `${item[1]}.h5ad`;
                    document.body.appendChild(a);
                    a.click();
                    
                    // 清理
                    window.URL.revokeObjectURL(url);
                    document.body.removeChild(a);
                  })
                  .catch(error=>{
                    console.error('Error fetching blob:', error);
                  })
                }}
              />
            </div>
          }>
          <List.Item.Meta
            title={
              <div>
                {item[0]}. <DatabaseOutlined /> {item[attr['title']]}
              </div>
            }
            description={
              <div>
                <b>Contributors:</b> {item[attr['contributors']]} <br />
                <TextCollapse
                  text={item[attr['summary']]}
                  header={<b>Summary: </b>}
                />
                <Space>
                {item[attr['pmid']] ? (
                  <div>
                    <b>PMID: </b> {item[attr['pmid']]}
                  </div>
                ) : (
                  ' '
                )}
                {
                  (  item[attr.spots] ? (
                    <div>
                      <Button type="primary" size="small" shape="round">
                        {AddCommas(item[attr.spots])} spots
                      </Button>
                    </div>
                  ) : (
                    ' '
                  ))
                }
                </Space>
              </div>
            }
          />
        </List.Item>
      )}
    />
  )
}
export default DatasetList
