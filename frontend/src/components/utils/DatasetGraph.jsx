import React, { useEffect, useRef, useState } from 'react';
import NeoVis from 'neovis.js';
import { Spin } from 'antd';
import loadingTips from '../charts/LoadingTip';
import Loading from '../charts/Loading';
const { enterLoading, quitLoading } = Loading


const DatasetGraph = ({ config }) => {
  const container = useRef(null);
  const [loadings, setLoadings] = useState([false]) // loadings for [network]
  const [currTip, setCurrTip] = useState(loadingTips[0])
  let vis = null
  let newConfig = {
    containerId: 'viz',
    neo4j:{
      serverUrl: "bolt://localhost:7687",
      serverUser: "neo4j",
      serverPassword: "biorzh123456",
    },
    labels: {
      ST: {
        label: 'id',  
        shape: 'circle',
        group: 'tissues',
      },
      SC: {
        label: 'id', 
        shape: 'text',
        group: 'tissues',
      },
    },
    relationships: {
      PAIR:{
          value: 5,  //String：线段粗细，用作边缘厚度的属性名。默认为1。
          caption: true,  //Boolean：如果设置为true，则关系类型将显示为边缘标题。或String：用作边缘标题的属性名。
          font: { size: 14, color: '#606266' }  // 关系节点文字大小颜色
      },
    },
    initialCypher: "MATCH p=()-[r:PAIR]->() RETURN p LIMIT 100"
  }

  useEffect(() => {
    enterLoading(0, setLoadings)
    setCurrTip(loadingTips[0])
    vis = new NeoVis(newConfig);
    vis.render()
    quitLoading(0, setLoadings)
  }, [config]);
 
  return (<Spin spinning={loadings[0]} size="large" tip={currTip}>
    <div ref={container} id='viz' style={{height:500}} />
    </Spin>
  );
};
 
export default DatasetGraph;