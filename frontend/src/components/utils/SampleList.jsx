import {
  LikeOutlined,
  MessageOutlined,
  StarOutlined,
  DatabaseOutlined,
  ExperimentOutlined,
  SwitcherOutlined,
  YuqueOutlined,
  FileTextOutlined,
  EyeOutlined,
} from '@ant-design/icons'
import React, { useRef, useState } from 'react'
import { Avatar, List, Space, Checkbox } from 'antd'

var data = [
  {
    key: 1,
    href: 'https://ant.design',
    title:
      'Spatial transcriptomics combined with single-cell RNA-sequencing unravels the complex inflammatory cell network in atopic dermatitis [ST]',
    avatar: `https://xsgames.co/randomusers/avatar.php?g=pixel&key=${0}`,
    description: (
      <div>
        <b>Contributors:</b> zhihan Ruan <br />
        Spatial Transcriptomics (Visium Spatial Gene Expression; 10X Genomics)
        of human skin biopsy samples. 14 skin biopsy samples were collected from
        7 atopic dermatitis patients, 1 from lesional and 1 from non-lesional
        skin per patient. In addition, 6 skin biopsy samples were collected from
        6 healthy donors to serve as a control.
      </div>
    ),
  },
]

var adata = Array.from({
  length: 23,
}).map((_, i) => ({
  key: i,
  href: 'https://ant.design',
  title: `ant design part ${i}`,
  avatar: `https://xsgames.co/randomusers/avatar.php?g=pixel&key=${i}`,
  description:
    'Ant Design, a design language for background applications, is refined by Ant UED Team.',
  content:
    'We supply a series of design principles, practical patterns and high quality design resources (Sketch and Axure), to help people create their product prototypes beautifully and efficiently.',
}))

data = [...data, ...adata]
const IconText = ({ icon, text }) => (
  <Space>
    {React.createElement(icon)}
    {text}
  </Space>
)

const SampleList = () => {
  const [selectedRowKeys, setSelectedRowKeys] = useState([])
  const rowSelection = {
    selectedRowKeys,
    onChange: (keys) => setSelectedRowKeys(keys),
  }
  return (
    <List
      itemLayout="vertical"
      size="large"
      rowSelection={rowSelection}
      pagination={{
        showSizeChanger: true,
        onChange: (page) => {
          console.log(page)
        },
        pageSize: 5,
      }}
      dataSource={data}
      renderItem={(item) => (
        <List.Item
          key={item.key}
          actions={[
            <IconText
              icon={YuqueOutlined}
              text="Homo Sapiens"
              key="list-vertical-star-o"
            />,
            <IconText
              icon={SwitcherOutlined}
              text="Skin"
              key="list-vertical-like-o"
            />,
            <IconText
              icon={ExperimentOutlined}
              text="10x Visium"
              key="list-vertical-message"
            />,
          ]}
          extra={
            <div>
              <br />
              <IconText
                icon={StarOutlined}
                text="Quality"
                key="list-vertical-star-o"
              />
              <br />
              <br />
              <IconText
                icon={FileTextOutlined}
                text="Description"
                key="list-vertical-star-o"
              />
              <br />
              <br />
              <IconText
                icon={EyeOutlined}
                text="Visualization"
                key="list-vertical-star-o"
              />
            </div>
          }>
          <List.Item.Meta
            title={
              <div>
                {item.key}. <DatabaseOutlined /> {item.title}
              </div>
            }
            description={item.description}
          />
          {item.content}
        </List.Item>
      )}
    />
  )
}
export default SampleList
