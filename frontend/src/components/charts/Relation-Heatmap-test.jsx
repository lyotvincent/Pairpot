import React, { useRef, useEffect, useState, useImperativeHandle } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import strokeColor from '../theme/strokeColor'
import { Space, Spin, Switch, Select, Progress } from 'antd'
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
import H5adLoader from '../utils/H5adLoader'
import Loading from './Loading'
import loadingTips from './LoadingTip'
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

const RelHeat = ({ spfile, scfile, setCompLoad, onRef, title, height, width, margin, progress }) => {
  const chartRef = useRef(null) // get current DOM container
  const commandRef = useRef('')
  const [action, setAction] = useState(0)
  const [Init, setInit] = useState(false)
  const [loadings, setLoadings] = useState([false]) // loadings for [heatmap]

  const [currentData, setCurrentData] = useState(null)
  const [currentInt, setCurrentInt] = useState(null)
  const [currentCell, setCurrentCell] = useState(null)
  const [currMethod, setCurrMethod] = useState(null)
  const [currCluster, setCurrCluster] = useState(null)

  const [methods, setMethods] = useState([])
  const [myCells, setmyCells] = useState([])
  const [DataArray, setDataArray] = useState([])
  const [IntArray, setIntArray] = useState([])
  const [CellArray, setCellArray] = useState([])
  const [dataKeys, setDataKeys] = useState([])

  const [scHeatmap, setScHeatmap] = useState({})
  const [spHeatmap, setSpHeatmap] = useState({})
  const [key, setKey] = useState(true)
  const [currTip, setCurrTip] = useState(loadingTips[0])
  const ScH5adLoader = (file) => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = (event) => {
        try {
          console.log("Call ScH5adLoader")
          let h5info = H5adLoader(file, event, ['cpdb'])
          resolve(h5info)
        } catch (error) {
          reject(error)
        }
      }
      reader.onloadend = () => {
        //console.log("Load heatmap data finished.")
      }
      reader.onerror = (error) => {
        reject(error)
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
    "Loading": setLoadings, // set Loading status
    "Tip": setCurrTip, // set Loading Tips
  }))

  useEffect(() => {
    if (Init) {
      var myChart = echarts.getInstanceByDom(chartRef.current)

      // console.log(currCluster)
      // console.log(myCells)

      if (commandRef.current === "Reload") {
        enterLoading(0, setLoadings)
        setCompLoad((compLoad) => {
          let newCompLoad = { ...compLoad }
          newCompLoad['CPDBHeatmap'] = true
          return newCompLoad
        })
        ScH5adLoader(scfile).then((h5info) => {
          setScHeatmap(h5info)
        })
        ScH5adLoader(spfile).then((h5info) => {
          console.log("All data reloaded in CPDB.")
          setSpHeatmap(h5info)
          setMethods(h5info.methods)
          setCurrMethod(h5info.methods[0])
          setCurrCluster(h5info.cellType[0])
          setDataArray(h5info.dataArray)
          setCellArray(h5info.cellArray)
          setIntArray(h5info.intArray)
          setmyCells(h5info.cellType)
          setDataKeys(h5info.dataKeys)
          setCurrentData(h5info.dataArray[h5info.cellType[0]])
          setCurrentInt(h5info.intArray[h5info.cellType[0]])
          setCurrentCell(h5info.cellArray[h5info.cellType[0]])
          toggleAnno("Upload")
        }).catch(error => {
          console.error('Error fetching blob in LayerView:', error)
        })
      }
      if (commandRef.current === "Upload") {
        enterLoading(0, setLoadings)
        const interactions = currentInt
        const celltypes = currentCell
        // console.log(currentInt)
        // console.log(currentCell)
        myChart.setOption({
          tooltip: {},
          title: [
            {
              text: title,
              left: '40%',
              textStyle: {
                fontSize: 24,
              },
            },
          ],
          grid: [
            {
              top: '5%',
              left: '20%',
              right: '15%',
              bottom: '25%',
              axisLine: {
                lineStyle: {
                  color: '#fff',
                },
              },
              axisPointer: {
                lineStyle: {
                  color: '#ffbd67',
                },
              },
            },
          ],
          yAxis: [
            {
              gridIndex: 0,
              name: 'L-R Pairs',
              nameLocation: 'end',
              nameGap: 10,
              nameRotate: 0,
              nameTextStyle: {
                fontWeight: 'bold',
                fontSize: 16,
                padding: 1,
              },
              type: 'category',
              data: interactions,
              splitArea: {
                show: true,
                interval: 4,
              },
              axisTick: {
                alignWithLabel: true,
              },
              axisLabel: {
                interval: 0,
                rotate: 0,
                margin: 10,
              },
            },
          ],
          xAxis: [
            {
              gridIndex: 0,
              name: 'Cell-type Pairs',
              nameLocation: 'end',
              nameGap: 10,
              nameRotate: 0,
              nameTextStyle: {
                fontWeight: 'bold',
                fontSize: 16,
                padding: 1,
              },
              type: 'category',
              data: celltypes,
              splitArea: {
                show: true,
              },
              axisTick: {
                alignWithLabel: true,
              },
              axisLabel: {
                interval: 0,
                fontSize: 14,
                rotate: -90,
                margin: 12,
              },
            },
          ],
          dataset: [{
            dimensions: dataKeys,
            source: currentData[0].map((_, i) => currentData.map(row => row[i])),
          },
          {
            // 这个 dataset 的 index 是 `1`。
            transform: {
              type: 'filter',
              config: { dimension: "significant", value: "yes" },
            },
            fromDatasetIndex: 0,
          },
          {
            // 这个 dataset 的 index 是 `2`。
            transform: {
              type: 'filter',
              config: { dimension: "significant", value: "nan" },
            },
            fromDatasetIndex: 0,
          }
          ],
          series: [
            {
              type: 'scatter',
              xAxisIndex: 0,
              yAxisIndex: 0,
              datasetIndex: 1,
              encode: {
                x: 'celltype_group',
                y: 'interaction_group',
                tooltip: [0, 1, 2, 3],
              },
              // borderColor: (val) =>{
              //   return val[8]==='yes' ? 'red' : '#fff'
              // },
              // borderColor :"red",
              itemStyle: {
                borderColor: 'red',
                borderWidth: 1,
              },
              //categories: celltypes,
              emphasis: {
                itemStyle: {
                  borderColor: 'red',
                  borderWidth: 1,
                },
              },
            },
            {
              type: 'scatter',
              xAxisIndex: 0,
              yAxisIndex: 0,
              datasetIndex: 2,
              encode: {
                x: 'celltype_group',
                y: 'interaction_group',
                tooltip: [0, 1, 2, 3],
              },
              // symbolSize: (val) => {
              //   return val[7] * 12 // 根据scaled_means
              // },
              // borderColor: (val) =>{
              //   return val[8]==='yes' ? 'red' : '#fff'
              // },
              // borderColor :"red",
              itemStyle: {
                borderColor: '#fff',
                borderWidth: 1,
              },
              //categories: celltypes,
              emphasis: {
                itemStyle: {
                  borderColor: 'red',
                  borderWidth: 1,
                },
              },
            },
          ],
          visualMap: [
            {
              left: 'right',
              top: '23%',
              orient: 'vertical',
              calculable: true,
              dimension: 'scaled_means',
              itemWidth: 20,
              itemHeight: 120,
              inRange: {
                // 设置映射的颜色范围为蓝色到红色之间
                symbolSize: [0, 18],
                color: ["#66CDAA", "#FFD700"],
              },
              text: ['scaled_means'],
              textGap: 20,
              textStyle: {
                fontSize: 14,
              },
              padding: 8,
              precision: 1,
              label: {
                formatter: '{value} %',
              },
              min: 0,
              max: 1,
              seriesIndex: [0, 1],
              handleIcon: 'circle',
              handleSize: '100%',
            },

          ],
        })
        quitLoading(0, setLoadings)
      }
      else if (commandRef.current === "Setting") {
        let interactions = currentInt
        let celltypes = currentCell
        myChart.setOption({
          dataset: [{
            dimensions: dataKeys,
            source: currentData[0].map((_, i) => currentData.map(row => row[i])),
          },
          {
            // datasetIndex is 1
            transform: {
              type: 'filter',
              config: { dimension: "significant", value: "yes" },
            },
            fromDatasetIndex: 0,
          },
          {
            // datasetIndex is 2
            transform: {
              type: 'filter',
              config: { dimension: "significant", value: "nan" },
            },
            fromDatasetIndex: 0,
          }

          ],
          xAxis: [
            {
              gridIndex: 0,
              data: celltypes,
            },
          ],
          yAxis: [
            {
              gridIndex: 0,
              data: interactions,
            }
          ],
        })
      }
    } else {
      var myChart = echarts.init(chartRef.current) //init the echart container
      enterLoading(0, setLoadings)
      setCurrTip(loadingTips[0])
      myChart.setOption({
        textStyle: {
          fontFamily: 'Arial',
        },
        tooltip: {},
        title: [
          {
            text: title,
            left: '40%',
            textStyle: {
              fontSize: 24,
            },
          },
        ],
        grid: [
          {
            top: '5%',
            left: '20%',
            right: '15%',
            bottom: '25%',
            axisLine: {
              lineStyle: {
                color: '#fff',
              },
            },
            axisPointer: {
              lineStyle: {
                color: '#ffbd67',
              },
            },
          },
        ],
        yAxis: [
          {
            gridIndex: 0,
            name: 'LR Pairs',
            nameLocation: 'end',
            nameGap: 70,
            nameRotate: 90,
            nameTextStyle: {
              fontWeight: 'bold',
              fontSize: 16,
              padding: 1,
            },
            type: 'category',
            splitArea: {
              show: true,
              interval: 4,
            },
            axisTick: {
              alignWithLabel: true,
            },
            axisLabel: {
              interval: 0,
              rotate: 0,
              margin: 10,
            },
          },
        ],
        xAxis: [
          {
            gridIndex: 0,
            name: 'Cell-type Pairs',
            nameLocation: 'middle',
            nameGap: 70,
            nameRotate: 0,
            nameTextStyle: {
              fontWeight: 'bold',
              fontSize: 16,
              padding: 1,
            },
            type: 'category',
            splitArea: {
              show: true,
            },
            axisTick: {
              alignWithLabel: true,
            },
            axisLabel: {
              interval: 0,
              fontSize: 14,
              rotate: -90,
              margin: 12,
            },
          },
        ],
        visualMap: [
          {
            left: 'right',
            top: '23%',
            orient: 'vertical',
            calculable: true,
            dimension: 'scaled_means',
            itemWidth: 20,
            itemHeight: 120,
            inRange: {
              // 设置映射的颜色范围为蓝色到红色之间
              symbolSize: [0, 18],
              color: ["#66CDAA", "#FFD700"],
            },
            text: ['scaled_means'],
            textGap: 20,
            textStyle: {
              fontSize: 14,
            },
            padding: 8,
            precision: 1,
            label: {
              formatter: '{value} %',
            },
            min: 0,
            max: 1,
            seriesIndex: [0, 1],
            handleIcon: 'circle',
            handleSize: '100%',
          },

        ],
      })
    }
    setInit(true)
  }, [action, title])


  // const buttonStyle = {
  //   backgroundColor: '#90EE90',
  //   color: 'white',
  //   width: '100px',
  //   height: '50px',
  //   border: 'none',
  //   borderRadius: '5px',
  //   margin: '5px'
  // }
  // ['Astrocytes', 'Interneurons', 'Neuroblasts1', 'Neuroblasts2', 'Neurons1', 'Neurons2',
  //   'Neurons3', 'Platelets', 'Pyramidal cells', 'Unknown']
  return (
    <Space>
      <Spin
        spinning={loadings[0]}
        size="large"
        tip={
          <div>{currTip}
            <Progress percent={progress} strokeColor={strokeColor} size={[300, 15]} />
          </div>
        }>
        <div
          ref={chartRef}
          className="chart"
          //the target DOM container needs height and width
          style={{ height: height, width: width, margin: margin }}></div>
      </Spin>
      <div style={{ marginLeft: -50 }}>
        <Space direction='vertical'>
          <Switch checkedChildren="Spatial"
            unCheckedChildren="Single-cell"
            value={key}
            onChange={(value) => {
              setKey(value)
              let h5info = value ? spHeatmap : scHeatmap
              setMethods(h5info.methods)
              setCurrMethod(h5info.methods[0])
              setCurrCluster(h5info.cellType[0])
              setDataArray(h5info.dataArray)
              setCellArray(h5info.cellArray)
              setIntArray(h5info.intArray)
              setmyCells(h5info.cellType)
              setDataKeys(h5info.dataKeys)
              setCurrentData(h5info.dataArray[h5info.cellType[0]])
              setCurrentInt(h5info.intArray[h5info.cellType[0]])
              setCurrentCell(h5info.cellArray[h5info.cellType[0]])
              toggleAnno("Upload")
            }} />
          <div>
            <span>Cluster: </span>
            <Select
              labelInValue
              placeholder={`(${myCells[0]})`}
              style={{
                width: '100%', margin: '2px'
              }}
              value={currCluster}
              options={myCells.map((cell, index) => ({ value: index, label: cell }))}
              onChange={(target) => {
                let selectedIndex = target.value
                setCurrentData(DataArray[myCells[selectedIndex]])
                setCurrentInt(IntArray[myCells[selectedIndex]])
                setCurrentCell(CellArray[myCells[selectedIndex]])
                setCurrCluster(myCells[selectedIndex])
                toggleAnno("Setting")
              }}
            />
          </div>
          <div>
            <span>Method: </span>
            <Select
              labelInValue
              placeholder={`(${methods[0]})`}
              style={{
                width: '100%', margin: '2px'
              }}
              value={currMethod}
              options={methods.map((cell, index) => ({ value: index, label: cell }))}
              onChange={(target) => {
                let selectedIndex = target.value
                let h5info = key ? spHeatmap : scHeatmap

                let prefix = ''
                if (methods[selectedIndex] !== 'CellphoneDB') {
                  prefix = `${methods[selectedIndex]}_`
                }
                setCurrMethod(methods[selectedIndex])
                setDataArray(h5info[`${prefix}dataArray`])
                setCellArray(h5info[`${prefix}cellArray`])
                setIntArray(h5info[`${prefix}intArray`])
                setmyCells(h5info[`${prefix}cellType`])
                setDataKeys(h5info[`dataKeys`])
                setCurrentData(h5info[`${prefix}dataArray`][h5info[`${prefix}cellType`][0]])
                setCurrentInt(h5info[`${prefix}intArray`][h5info[`${prefix}cellType`][0]])
                setCurrentCell(h5info[`${prefix}cellArray`][h5info[`${prefix}cellType`][0]])
                setCurrCluster(h5info[`${prefix}cellType`][0])
                toggleAnno("Upload")
              }}
            />
          </div>
        </Space>
      </div>
    </Space>
  )
}

RelHeat.defaultProps = {
  title: 'L-R Pairs Heatmap',
  height: '60rem',
  width: '60rem',
  margin: '2rem',
}

RelHeat.propTypes = {
  spfile: PropTypes.object,
  scfile: PropTypes.object,
  setCompLoad: PropTypes.func,
  title: PropTypes.string,
  onRef: PropTypes.any,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
  progress: PropTypes.number,
}

export default RelHeat
