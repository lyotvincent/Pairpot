import React, { useRef, useEffect, useState, useImperativeHandle } from 'react'
import * as echarts from 'echarts/core'
import PropTypes from 'prop-types'
import H5adLoader from '../utils/H5adLoader'
import {
  GraphicComponent,
  GridComponent,
  LegendComponent,
  TooltipComponent,
  VisualMapComponent,
} from 'echarts/components'
import { LineChart, LinesChart } from 'echarts/charts'
import { UniversalTransition } from 'echarts/features'
import { CanvasRenderer } from 'echarts/renderers'
import Loading from './Loading'
import { Spin, Switch } from 'antd'
echarts.use([
  GraphicComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  LineChart,
  CanvasRenderer,
  UniversalTransition,
  VisualMapComponent,
  LinesChart,
])

const { enterLoading, quitLoading } = Loading



const NetworkRelation = ({ spfile, scfile, onRef, title, height, width, margin }) => {
  const chartRef = useRef(null)
  const commandRef = useRef('')
  const [action, setAction] = useState(0)
  const [key, setKey] = useState(true)
  const [Init, setInit] = useState(false)
  const [loadings, setLoadings] = useState([false]) // loadings for [network]
  const [graph, setGraph] = useState(null)
  const [scgraph, setSCGraph] = useState(null)
  const [spgraph, setSPGraph] = useState(null)

  const ScH5adLoader = (file) => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = (event) => {
        try {
          let h5info = H5adLoader(file, event, ['net'])
          h5info = {
            links: h5info.links,
            categories: h5info.categories,
            nodes: h5info.nodes
          }
          h5info['links'].forEach(function (link) {
            link['lineStyle'] = {
              width: link.value / 6,
              curveness: link.source === link.target ? 0.9 : 0, // 自连边设置为0.8，其他边设置为0.3
              label: {
                show: link.source === link.target ? false : true
              }
            }
          })
          resolve({
            links: h5info.links,
            categories: h5info.categories,
            nodes: h5info.nodes
          })
        } catch (error) {
          reject(error)
        }
      }
      reader.onloadend = () => {
        //console.log("Load network data finished.")
      }
      reader.onerror = (error) => {
      }
      reader.readAsArrayBuffer(file)
    })

  }

  const toggleAnno = (command) => {
    commandRef.current = command
    setAction(action + 1)
  }

  useImperativeHandle(onRef, () => ({  // explode trigger for parent components
    "Trigger": toggleAnno, // Trigger for useEffect
  }))

  useEffect(() => {
    if (Init) {
      var myChart = echarts.getInstanceByDom(chartRef.current)
      if (commandRef.current === "Reload") {
        enterLoading(0, setLoadings)
        scfile.then((loadedFile) => {
          ScH5adLoader(loadedFile).then((res) => {
            setSCGraph(res)
          })
        })
        spfile.then((loadedFile) => {
          ScH5adLoader(loadedFile).then((res) => {
            console.log("All data reloaded in Network.")
            setSPGraph(res)
            setGraph(res)
            toggleAnno("Upload")
          })
        })
      }
      if (commandRef.current === "Upload") {
        enterLoading(0, setLoadings)
        myChart.setOption({
          legend: [
            {
              data: graph.categories.map((a) => {
                return a.name
              })
            }
          ],
          series: [
            {
              name: 'interaction-count',
              type: 'graph',
              layout: 'force',
              data: graph.nodes,
              links: graph.links,
              categories: graph.categories,
              roam: false,
              focusNodeAdjacency: true,
              // 力引导布局相关的配置项
              force: {
                gravity: 1,
                repulsion: 100,
                edgeLength: [500, 200]
              },
              label: {
                show: true
              },
              edgeSymbol: ['circle', 'arrow'],
              edgeSymbolSize: [4, 10],
              edgeLabel: {
                fontSize: 20
              },
              label: {
                show: true,
                position: 'midlle',
                formatter: '{b}'
              },
              labelLayout: {
                hideOverlap: true
              },
              scaleLimit: {
                min: 0.4,
                max: 2
              },
              lineStyle: {
                color: '#B0E0E6',
                curveness: 0,
                opacity: 1
                //   width: 
              },
            }
          ],
        }
        )
        quitLoading(0, setLoadings)
      }
    } else {
      var myChart = echarts.init(chartRef.current) //init the echart container
      enterLoading(0, setLoadings)
      scfile.then((loadedFile) => {
        ScH5adLoader(loadedFile).then((res) => {
          setSCGraph(res)
        })
      })
      spfile.then((loadedFile) => {
        ScH5adLoader(loadedFile).then((res) => {
          console.log("All data loaded in Network.")
          setSPGraph(res)
          setGraph(res)
          toggleAnno("Upload")
        })
      })
      myChart.setOption({
        textStyle: {
          fontFamily: 'Arial',
        },
        tooltip: {},
        legend: {
          animation: false,
          orient: 'horizontal',
          top: '5%',
          padding: 10,
          textStyle: {
            fontSize: 16,
          },
        },
        title: [
          {
            text: title,
            left: 'center',
            textStyle: {
              fontSize: 24,
            },
          },
        ],
        animation: false,
        grid: [{
          top: '0%',
          bottom: '0%',
        }],
      })
      setInit(true)
    }
  }, [action])


  return (
    <Spin spinning={loadings[0]} size="large" tip="Loading">
      <Switch checkedChildren="Spatial"
        unCheckedChildren="Single-cell"
        value={key}
        onChange={(value) => {
          setKey(value)
          setGraph(value ? spgraph : scgraph)
          toggleAnno("Upload")
        }} />
      <div
        ref={chartRef}
        className="chart"
        //the target DOM container needs height and width
        style={{ height: height, width: width, margin: margin }}>
      </div>
    </Spin>
  )

  // return <div ref={chartRef} style={{ width: '100%', height: '400px' }} />;

}

// some settings
NetworkRelation.defaultProps = {
  title: 'Cell-Cell Interaction Network',
  height: '45rem',
  width: '40rem',
  margin: '0rem',
}

NetworkRelation.propTypes = {
  spfile: PropTypes.object,
  scfile: PropTypes.object,
  onRef: PropTypes.any,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default NetworkRelation

