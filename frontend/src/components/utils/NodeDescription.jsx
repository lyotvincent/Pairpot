import React, { useRef, useState, useEffect } from 'react'
import { Descriptions, Modal, Space, Tooltip,Button,Row,Col } from 'antd'
import { FileTextOutlined } from '@ant-design/icons'
import { useNavigate } from 'react-router-dom'
import axios from 'axios'
import { useQuery } from 'react-query'



const DatasetDescription = ({ descCol, descInfo, text, placement }) => {
  const [modal1Open, setModal1Open] = useState(false)
  const [descItem, setDescItem] = useState([])
  const navigate = useNavigate()

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

  const aItems = ['accessions']

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

    const OriginResponse = useQuery({
      queryKey: ['db'],
      queryFn: () => axios.get('/api/datasets').then((response) => {
        return response.data
      }).catch((error) => {
        console.log(error)
      }),
      staleTime: Infinity,
      retry: false,
      refetchOnWindowFocus: false,
    })

    const GetscItem = (scid) => {
      return OriginResponse.data?.data.find((item) => item[1] === scid)
    }

  useEffect(() => {
    let items = []
    if (descInfo.length > 0) {
      for (var i = 0; i < descCol.length; i++) {
        let useSpan = midItemSpan
        useSpan = lgItems.includes(descCol[i]) ? lgItemSpan : useSpan
        useSpan = smItems.includes(descCol[i]) ? smItemSpan : useSpan
        let useInfo = aItems.includes(descCol[i]) ? (
          <a href={descInfo[i]}>{descInfo[i]}</a>
        ) : (
          descInfo[i]
        )
        if (!excludeItems.includes(descCol[i])) {
          items.push({
            label: descCol[i],
            children: useInfo,
            span: useSpan,
          })
        }
      }
      setDescItem(items)
    }
  }, [descCol, descInfo])
  return (
    <div>
      {/* <Tooltip
        placement={placement}
        title="Show Descriptions"
        arrow={false}
        align={'center'}>
        <Space
          onClick={() => setModal1Open(true)}
          style={{ cursor: 'pointer' }}>
          <FileTextOutlined />
          {text}
        </Space>
      </Tooltip>
      <Modal
        title="Dataset Descriptions"
        open={modal1Open}
        onOk={() => setModal1Open(false)}
        onCancel={() => setModal1Open(false)}
        style={{ top: '15%' }}
        width={1000}> */}
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
      <Row justify="center" style={{marginTop:'20px'}}>
      <Col>
      <Button 
        type='primary'
        onClick={() => {
          let item = descInfo
          let dataCol = descCol
          console.log(dataCol)
          let values = Object.fromEntries(
              dataCol.map((k, i) => [k, item[i]])
          )
          let scid = values['has_paired']
          let scitem = GetscItem(scid)
          let state = {
          st: values
          }
          if (typeof scitem !== 'undefined') {
          let scvalues = Object.fromEntries(
              dataCol.map((k, i) => [k, scitem[i]])
          )
          state['sc'] = scvalues
          }
          else{
          console.log("No paired sc data found.")
          }
          navigate('/browse', { state: state })
      }}>
        visualization</Button>
        </Col>
        </Row>
      {/* </Modal> */}
    </div>
  )
}
export default DatasetDescription