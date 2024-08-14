import React, { useEffect, useRef, useState } from 'react';
import { NeoVis } from 'neovis.js';
import { Spin, Card, Col, Button, Row, Space } from 'antd';
import loadingTips from '../charts/LoadingTip';
import Loading from '../charts/Loading';
const { enterLoading, quitLoading } = Loading


const DatasetGraph = ({ config }) => {
  const container = useRef(null);
  const [loadings, setLoadings] = useState([false]) // loadings for [network]
  const [currTip, setCurrTip] = useState(loadingTips[0])
  const [rerender, setRerender] = useState(false)
  const visRef = useRef(null)
  let newConfig = {
    containerId: 'viz',
    neo4j: {
      serverUrl: "bolt://localhost:7687",
      serverUser: "neo4j",
      serverPassword: "biorzh123456",
    },
    labels: {
      ST: {
        label: 'id',
        font: {
          size: 18,
        },
        group: 'tissues',
      },
      SC: {
        label: 'id',
        group: 'tissues',
      },
    },
    hierarchical: true, // 节点显示方式 是否启用分层布局后
    // 分层结构或者默认 "hubsize"（默认）和"directed".
    hierarchical_sort_method: 'hubsize',
    visConfig: {
      nodes: {
        shape: 'dot',
        scaling: {
          max: 10,
          min: 10
        },
      },
      edges: {
        arrows: {
          to: { enabled: true }
        }
      },
    },
    relationships: {
      PAIR: {
        value: 5,  //String：线段粗细，用作边缘厚度的属性名。默认为1。
        caption: true,  //Boolean：如果设置为true，则关系类型将显示为边缘标题。或String：用作边缘标题的属性名。
        font: { size: 14, color: '#606266' }  // 关系节点文字大小颜色
      },
    },
  }
  useEffect(() => {
    enterLoading(0, setLoadings)
    setCurrTip(loadingTips[0])
    visRef.current = new NeoVis(newConfig)
    visRef.current.renderWithCypher("MATCH p=()-[r:PAIR]->() RETURN p LIMIT 25")
    quitLoading(0, setLoadings)
  }, []);

  useEffect(() => {
    if (config) {
      
      visRef.current?.renderWithCypher(`MATCH p=(n:ST ${config})-[r:PAIR]->() RETURN p LIMIT 100`)
    }
  }, [config])


  return (<Spin spinning={loadings[0]} size="large" tip={currTip}>
    <Row>
      <Col span={16} offset={1}>
        <Card>
          <div ref={container} id='viz' style={{ width: "100%", height: 500 }} />
        </Card>
      </Col>
      <Col span={6} offset={1}>
        <Space direction='vertical'>
          <Button type='primary' onClick={() => {
            visRef.current?.renderWithCypher("MATCH p=()-[r:PAIR]->() RETURN p LIMIT 100")
          }}> Display All
          </Button>
          <Button type='primary' onClick={() => {
            visRef.current?.renderWithCypher("MATCH p=(n:ST)-[r:PAIR]->(m:SC) RETURN p LIMIT 100")
          }}>
            Display SC-SP pairs
          </Button>
          <Button type='primary' onClick={() => {
            visRef.current?.renderWithCypher("MATCH p=(n:ST)-[r:PAIR]->(m:ST) RETURN p LIMIT 100")
          }}> 
          Display SRT Rings 
          </Button>
          <Button type='primary' onClick={() => {
            visRef.current?.stabilize()
          }}> 
          Stabilize
          </Button>
        </Space>
      </Col>
    </Row>
  </Spin>
  );
};

export default DatasetGraph;