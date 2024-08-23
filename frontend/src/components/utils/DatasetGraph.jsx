import React, { useEffect, useRef, useState } from 'react';
import { NeoVis } from 'neovis.js';
import { Spin, Card, Col, Button, Row, Space, Modal, Descriptions } from 'antd';
import loadingTips from '../charts/LoadingTip';
import Loading from '../charts/Loading';
import ToggleAccordion from './ToggleAccordion';
import { ChevronDoubleLeft } from 'react-bootstrap-icons';
import DatasetDescription from './NodeDescription'
const { enterLoading, quitLoading } = Loading


const DatasetGraph = ({ config, src, col }) => {
  const container = useRef(null);
  const [loadings, setLoadings] = useState([false]) // loadings for [network]
  const [currTip, setCurrTip] = useState(loadingTips[0])
  // const [rerender, setRerender] = useState(false)
  const visRef = useRef(null)

  const [modalVisible, setModalVisible] = useState(false); // 控制模态框的显示
  const [nodeData, setNodeData] = useState({}); // 保存点击节点的数据

  let newConfig = {
    containerId: 'viz',
    neo4j: {
      serverUrl: "bolt://123.207.205.5:17687",
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
    
    // 监听节点点击事件
    visRef.current.registerOnEvent('clickNode', (event) => {
      // const clickedNodeData = event.node.raw.properties; // 获取节点的属性数据
      const id = event.node.label
      // 在src中找到完整的id的数据
      const clickedNodeData = src.find(s => Array.isArray(s) && s[1] === id)
      setNodeData(clickedNodeData); // 设置节点数据
      setModalVisible(true); // 打开模态框
    });

    quitLoading(0, setLoadings)
  }, []);

  useEffect(() => {
    if (config) {

      visRef.current?.renderWithCypher(`MATCH p=(n:ST ${config})-[r:PAIR]->() RETURN p LIMIT 100`)
    }
  }, [config])


  return (<Spin spinning={loadings[0]} size="large" tip={currTip}>
    <ToggleAccordion header={
      <h3 id="Association Networks">{"Association Networks"}
      </h3>}>
      {<div>
        Pairpot provides association networks to establish connections among different single-cell and SRT studies。
        In the association networks, each node denotes a dataset and its color denotes the dataset's species.
        Nodes with self-loops in the network represent datasets that contain SC-SP pairs.
        Paired single-cell and SRT datasets are linked by undirected edges.
        <br/>
        Users can click the nodes to explore detailed information and access the analysis modules.
        By default, only 25 pairs would be displayed initially for the system fluency, if users would like to browse all the pairs, please click the button `Display All`.
      </div>}
    </ToggleAccordion>
    <br />
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

    {/* 节点信息的模态框 */}
    <Modal
        title="Node Details"
        open={modalVisible}
        onCancel={() => setModalVisible(false)}
        footer={null}
        width={1200}
        // height={400}
      >
        {/* <Descriptions bordered column={2}>
          {Object.entries(nodeData)
            .filter(([key]) => ['id', 'species', 'technologies', 'tissues', 'title'].includes(key))  // 过滤特定字段
            .map(([key, value]) => (
              <Descriptions.Item label={key} key={key}>
                {value}
              </Descriptions.Item>
            ))}
        </Descriptions> */}
        <DatasetDescription
          descCol = {col}
          descInfo = {nodeData}
          text={'Descriptions'}
          placement={'left'}
        >
        </DatasetDescription>
      </Modal>
  </Spin>
  );
};

export default DatasetGraph;