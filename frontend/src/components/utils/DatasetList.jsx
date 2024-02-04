import {
  DatabaseOutlined,
  ExperimentOutlined,
  SwitcherOutlined,
  YuqueOutlined,
  EyeOutlined,
  LinkOutlined,
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
//import { response } from 'express'

const IconText = ({ icon, text, attr, placement = 'bottom' }) => (
  <Space>
    <Tooltip placement={placement} title={attr} arrow={false} align={'center'}>
      {React.createElement(icon)}
    </Tooltip>
    {text}
  </Space>
)

const AddCommas = (num) => {
  return num.toLocaleString().replace(/^0+/, '')
}

const DatasetList = ({ src, col }) => {
  const [dataSrc, setDataSrc] = useState([])
  const [attr, setAttr] = useState([])
  const [dataCol, setDataCol] = useState({})
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
            />,
            <Rate defaultValue={2} />,
          ]}
          extra={
            <div>
              {item[attr.spots] ? (
                <div>
                  <br />
                  <Button type="primary" size="small" shape="round">
                    {AddCommas(item[attr.spots])} spots
                  </Button>
                </div>
              ) : (
                ' '
              )}
              <br />
              <DatasetDescription descCol={dataCol} descInfo={item} />
              <br />

              <IconText
                icon={EyeOutlined}
                text="Visualization"
                key="list-vertical-visualization-o"
                attr="Visualize this Dataset"
                placement="left"
              />
              <br />
              <br />
              <IconText
                icon={LinkOutlined}
                text="Paired Datasets"
                key="list-vertical-Paired-Datasets-o"
                attr="Show paired datasets"
                placement="left"
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
                {item[attr['pmid']] ? (
                  <div>
                    <b>PMID: </b> {item[attr['pmid']]}
                  </div>
                ) : (
                  ' '
                )}
              </div>
            }
          />
        </List.Item>
      )}
    />
  )
}
export default DatasetList
