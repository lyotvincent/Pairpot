import React, { useRef, useState, useEffect } from 'react'
import { Descriptions, Modal, Space, Tooltip } from 'antd'
import { FileTextOutlined } from '@ant-design/icons'
const items = [
  {
    label: 'Product',
    children: 'Cloud Database',
  },
  {
    label: 'Billing',
    children: 'Prepaid',
  },
  {
    label: 'Time',
    children: '18:00:00',
  },
  {
    label: 'Amount',
    children: '$80.00',
  },
  {
    label: 'Discount',
    span: {
      xl: 2,
      xxl: 2,
    },
    children: '$20.00',
  },
  {
    label: 'Official',
    span: {
      xl: 2,
      xxl: 2,
    },
    children: '$60.00',
  },
  {
    label: 'Config Info',
    span: {
      xs: 1,
      sm: 2,
      md: 3,
      lg: 3,
      xl: 2,
      xxl: 2,
    },
    children: (
      <>
        Data disk type: MongoDB
        <br />
        Database version: 3.4
        <br />
        Package: dds.mongo.mid
      </>
    ),
  },
  {
    label: 'Hardware Info',
    span: {
      xs: 1,
      sm: 2,
      md: 3,
      lg: 3,
      xl: 2,
      xxl: 2,
    },
    children: (
      <>
        CPU: 6 Core 3.5 GHz
        <br />
        Storage space: 10 GB
        <br />
        Replication factor: 3
        <br />
        Region: East China 1
      </>
    ),
  },
]
const DatasetDescription = ({ descCol, descInfo }) => {
  const [modal1Open, setModal1Open] = useState(false)
  const [descItem, setDescItem] = useState({})

  const lgItemSpan = {
    xl: 6,
    xxl: 6,
  }

  const midItemSpan = {
    xs: 1,
    sm: 2,
    md: 3,
    lg: 3,
    xl: 3,
    xxl: 3,
  }

  const smItemSpan = {
    xs: 1,
    sm: 1,
    md: 2,
    lg: 2,
    xl: 2,
    xxl: 2,
  }

  const lgItems = [
    'summary',
    'overall_design',
    'citation',
    'accessions',
    'title',
    'contributors',
    'contacts',
    'development_stages',
  ]
  const smItems = [
    'species',
    'tissues',
    'organ_parts',
    'cells',
    'spots',
    'genes',
    'sex',
    'technologies',
    '',
  ]

  const excludeItems = ['cell_types']

  useEffect(() => {
    var items = []
    if (descInfo.length > 0) {
      for (var i = 0; i < descCol.length; i++) {
        var useSpan = midItemSpan
        useSpan = lgItems.includes(descCol[i]) ? lgItemSpan : useSpan
        useSpan = smItems.includes(descCol[i]) ? smItemSpan : useSpan
        if (!excludeItems.includes(descCol[i])) {
          items.push({
            label: descCol[i],
            children: descInfo[i],
            span: useSpan,
          })
        }
      }
      setDescItem(items)
    }
  }, [descCol, descInfo])
  return (
    <div>
      <Tooltip
        placement="left"
        title="Show Descriptions"
        arrow={false}
        align={'center'}>
        <Space
          onClick={() => setModal1Open(true)}
          style={{ cursor: 'pointer' }}>
          <FileTextOutlined /> Descriptions
        </Space>
      </Tooltip>
      <Modal
        title="Dataset Descriptions"
        open={modal1Open}
        onOk={() => setModal1Open(false)}
        onCancel={() => setModal1Open(false)}
        style={{ top: '15%' }}
        width={1000}>
        <Descriptions
          bordered
          size="small"
          column={{
            xs: 2,
            sm: 4,
            md: 6,
            lg: 6,
            xl: 6,
            xxl: 6,
          }}
          items={descItem}
        />
      </Modal>
    </div>
  )
}
export default DatasetDescription
