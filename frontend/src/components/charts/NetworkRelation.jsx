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
import strokeColor from '../theme/strokeColor'
import { LineChart, LinesChart } from 'echarts/charts'
import { UniversalTransition } from 'echarts/features'
import { CanvasRenderer } from 'echarts/renderers'
import Loading from './Loading'
import { Col, Spin, Switch, Select, Button, Progress, Space } from 'antd'
import loadingTips from './LoadingTip'
import Axis from './Axis'
import { color } from 'echarts'
import vega_20 from '../theme/vega_20'
import { useFetcher } from 'react-router-dom'
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



const NetworkRelation = ({ spfile, scfile, setCompLoad, onRef, title, height, width, margin, progress }) => {
  const chartRef = useRef(null)
  const chartRef1 = useRef(null) // layer view + lasso view
  const commandRef = useRef('')
  const [action, setAction] = useState(0)
  const [key, setKey] = useState(true)
  const [Init, setInit] = useState(false)
  const [loadings, setLoadings] = useState([false]) // loadings for [network]
  const [graph, setGraph] = useState(null)
  const [scgraph, setSCGraph] = useState(null)
  const [spgraph, setSPGraph] = useState(null)
  const [currTip, setCurrTip] = useState(loadingTips[0])
  const [_data, _setData] = useState([])
  const [_data1, _setData1] = useState([])
  const [clusterOps, setClusterOps] = useState([])
  const [cellNum, setCellNum] = useState(0)
  const [clusterCur, setClusterCur] = useState({})
  const prevCluster = useRef('')
  const itemGroupRef = useRef([])
  const [batchName, setBatchName] = useState("batch")
  const [batchCur, setBatchCur] = useState({})
  const [batchOps, setBatchOps] = useState([])
  const [embedCur, setEmbedCur] = useState('')
  const [embedCur1, setEmbedCur1] = useState('')
  const [embedOps, setEmbedOps] = useState([])
  const [xInv, setxInv] = useState(false)
  const [yInv, setyInv] = useState(false)
  const symbolSizeRef = useRef('')
  const [itemSize, setItemSize] = useState(2)
  const symbolSizeRef1 = useRef('')
  const [itemSize1, setItemSize1] = useState(2)
  const [cellArray, setCellArray] = useState([])
  const [lfc, setLFC] = useState([])  // logfoldchanges
  const [mn, setMN] = useState([])  // marker names
  const [pa, setPA] = useState([]) // p-value_adj
  const [me, setME] = useState([]) // mean expression
  const [fr, setFR] = useState([]) // fracments
  const [den, setDen] = useState({}) // dendrogram
  const [arrLen, setArrLen] = useState(0) // Length of 1D lfc, mn and pval
  const [cluLen, setCluLen] = useState(0) // Length of clusters, arrLen / genes
  const geneNum = 2000
  const [lfc1, setLFC1] = useState([])  // logfoldchanges
  const [mn1, setMN1] = useState([])  // marker names
  const [pa1, setPA1] = useState([]) // p-value_adj
  const [me1, setME1] = useState([]) // mean expression
  const [fr1, setFR1] = useState([]) // fracments
  const [cluLen1, setCluLen1] = useState(0) // Length of clusters, arrLen / genes




  const setItemGroup = (source, annoIdx, type = 'categories') => {
    let itemGroup = []
    let annotations = []
    if (type === 'categories') {
      annotations = [
        ...new Set(
          source.map((item) => {
            return item[annoIdx]
          })
        ),
      ]

      annotations = annotations.sort()
      for (let anno of annotations) {
        itemGroup.push(
          source
            .filter((item) => {
              return item[annoIdx] === anno
            })
            .map((item) => item[item.length - 1])
        )
      }
    } else {
      itemGroup.push(source.map((item) => item[item.length - 1]))
    }

    itemGroupRef.current = itemGroup
    return annotations
  }


  const setAxis = (source, dims, xName0, yName0) => {
    // the axis of left chart
    let xIdx0 = dims.indexOf(xName0)
    let yIdx0 = dims.indexOf(yName0)
    let embd0_x = source.map((item) => {
      return item[xIdx0]
    })
    let embd0_y = source.map((item) => {
      return item[yIdx0]
    })

    return {
      xAxis: [
        {
          id: 0,
          gridIndex: 0,
          name: xName0,
          nameLocation: 'middle',
          nameGap: 23,
          nameTextStyle: {
            fontSize: 16,
          },
          axisLine: {
            onZero: false,
          },
          axisLabel: {
            fontSize: 14,
          },
          position: 'bottom',
          inverse: xInv,
          min: Math.ceil(Math.min.apply(null, embd0_x) - 1.5),
          max: Math.ceil(Math.max.apply(null, embd0_x) + 1.5),
        },
      ],
      yAxis: [
        {
          id: 0,
          gridIndex: 0,
          name: yName0,
          nameLocation: yInv ? 'start' : 'end',
          nameTextStyle: {
            fontSize: 16,
          },
          axisLine: {
            onZero: false,
          },
          axisLabel: {
            fontSize: 14,
          },
          position: 'left',
          inverse: yInv,
          min: Math.ceil(Math.min.apply(null, embd0_y) - 1.5),
          max: Math.ceil(Math.max.apply(null, embd0_y) + 1.5),
        },
      ],
      // xAxis3D: [{
      //   grid3DIndex: 0,
      //   name: xName0,
      //   nameLocation: 'middle',
      //   min: Math.ceil(Math.min.apply(null, embd0_x) - 1.5),
      //   max: Math.ceil(Math.max.apply(null, embd0_x) + 1.5),
      // }
      // ],
      // yAxis3D: [{
      //   grid3DIndex: 0,
      //   name: yName0,
      //   nameLocation: 'middle',
      //   min: Math.ceil(Math.min.apply(null, embd0_y) - 1.5),
      //   max: Math.ceil(Math.max.apply(null, embd0_y) + 1.5),
      // }
      // ],
    }
  }


  const setBatchDataset = (source, dims, annoName, annotations, batch, batchName) => {
    let _datasets = []
    _datasets.push({
      dimensions: dims,
      source: source,
    })
    if (batchName !== null) {
      if (annotations.length === 0) {
        // if no annotations, all the data are annotated by a label.
        _datasets.push({
          // 这个 dataset 的 index 是 `1`。

          transform: [
            {
              type: 'filter', //batchName is considered a column name in.obs
              config: { dimension: batchName, value: batch },
            },
            {
              type: 'sort',
              config: { dimension: annoName, order: 'desc' },
            },]
        })
      } else {
        // add datasets with batch filter
        for (let anno of annotations) {
          _datasets.push({
            transform: [
              {
                type: 'filter', //batchName is considered a column name in.obs
                config: { dimension: batchName, value: batch },
              },
              {
                type: 'filter',
                config: { dimension: annoName, value: anno },
              },]
          })
        }
        // add datasets without batch filter
        for (let anno of annotations) {
          _datasets.push({
            transform: [
              {
                type: 'filter',
                config: { dimension: annoName, value: anno },
              },]
          })
        }
      }
    } else {
      // add datasets without batch filter
      for (let anno of annotations) {
        _datasets.push({
          transform: [
            {
              type: 'filter',
              config: { dimension: annoName, value: anno },
            },]
        })
      }
      for (let anno of annotations) {
        _datasets.push({
          transform: [
            {
              type: 'filter',
              config: { dimension: annoName, value: anno },
            },]
        })
      }
    }
    return _datasets
  }

  const ScH5adLoader = (file, sp) => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = (event) => {
        try {
          let h5info = H5adLoader(file, event, ['net'])
          let h5info1 = H5adLoader(file, event, ['meta'])
          let h5info2 = H5adLoader(file, event, ['marker', 'den'])
          console.log(h5info2)
          if (sp) {
            _setData(h5info1.data)
            // top markers
            setMN(h5info2.names)
            setLFC(h5info2.lfc)
            setPA(h5info2.pval)
            setME(h5info2.expr)
            setFR(h5info2.frac)
            setDen({
              'icoord': h5info2.icoord,
              'dcoord': h5info2.dcoord,
              'ivl': h5info2.ivl,
              'leaves': h5info2.leaves,
              'color_list': h5info2.color_list,
              'leaves_color_list': h5info2.leaves_color_list
            })
            setArrLen(h5info2.names.length)
            setCluLen(h5info2.names.length / geneNum)            
          }
          else {
            _setData1(h5info1.data)
            setMN1(h5info2.names)
            setLFC1(h5info2.lfc)
            setPA1(h5info2.pval)
            setME1(h5info2.expr)
            setFR1(h5info2.frac)
            setCluLen1(h5info2.names.length / geneNum)
          }
          setClusterOps(h5info1.clusters)
          if (h5info1.clusters.map((item) => item.label).includes('batch'))
            setBatchName('batch')
          setEmbedOps(h5info1.embdOps)
          h5info = {
            links: h5info.links,
            categories: h5info.categories,
            nodes: h5info.nodes
          }
          if (sp) {
            // sort
            h5info.nodes.sort((a, b) => Number(a.category) - Number(b.category))
            h5info.categories.sort((a, b) => Number(a.name) - Number(b.name))
          }

          console.log(h5info)

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

  const setDataset = (source, dims, annoName, annotations) => {
    let _datasets = []
    _datasets.push({
      dimensions: dims,
      source: source,
    })
    if (annotations.length === 0) {
      // if no annotations, all the data are annotated by a label.
      _datasets.push({
        transform: {
          type: 'sort',
          config: { dimension: annoName, order: 'desc' },
        },
      })
    } else {
      for (let anno of annotations) {
        _datasets.push({
          transform: {
            type: 'filter',
            config: { dimension: annoName, value: anno },
          },
        })
      }
    }
    return _datasets
  }


  const setAxis1 = (
    source,
    dims,
    xName0,
    yName0,
    xName1 = null,
    yName1 = null
  ) => {
    // the axis of left chart
    let xIdx0 = dims.indexOf(xName0)
    let yIdx0 = dims.indexOf(yName0)
    let embd0_x = source.map((item) => {
      return item[xIdx0]
    })
    let embd0_y = source.map((item) => {
      return item[yIdx0]
    })

    // the axis of right chart. By default it's as same as the axis in left chart.
    let xIdx1 = xIdx0
    let yIdx1 = yIdx0
    let embd1_x = embd0_x
    let embd1_y = embd0_y
    if (xName1 !== null && yName1 !== null) {
      xIdx1 = dims.indexOf(xName1)
      yIdx1 = dims.indexOf(yName1)
      embd1_x = source.map((item) => {
        return item[xIdx1]
      })
      embd1_y = source.map((item) => {
        return item[yIdx1]
      })
    } else {
      xName1 = xName0
      yName1 = yName0
    }

    return {
      xAxis: [
        {
          id: 0,
          gridIndex: 0,
          name: xName0,
          nameLocation: 'middle',
          nameGap: 23,
          nameTextStyle: {
            fontSize: 16,
          },
          axisLine: {
            onZero: false,
          },
          axisLabel: {
            fontSize: 14,
          },
          position: 'bottom',
          inverse: xInv,
          min: Math.ceil(Math.min.apply(null, embd0_x) - 1.5),
          max: Math.ceil(Math.max.apply(null, embd0_x) + 1.5),
        },
        // {
        //   id: 1,
        //   gridIndex: 1,
        //   name: xName1,
        //   nameLocation: 'middle',
        //   nameGap: 23,
        //   nameTextStyle: {
        //     fontSize: 16,
        //   },
        //   axisLine: {
        //     onZero: false,
        //   },
        //   axisLabel: {
        //     fontSize: 14,
        //   },
        //   position: 'bottom',
        //   inverse: xInv,
        //   min: Math.ceil(Math.min.apply(null, embd1_x) - 1.5),
        //   max: Math.ceil(Math.max.apply(null, embd1_x) + 1.5),
        // },
      ],
      yAxis: [
        {
          id: 0,
          gridIndex: 0,
          name: yName0,
          nameLocation: yInv ? 'start' : 'end',
          nameTextStyle: {
            fontSize: 16,
          },
          axisLine: {
            onZero: false,
          },
          axisLabel: {
            fontSize: 14,
          },
          position: 'left',
          inverse: yInv,
          min: Math.ceil(Math.min.apply(null, embd0_y) - 1.5),
          max: Math.ceil(Math.max.apply(null, embd0_y) + 1.5),
        },
        // {
        //   id: 1,
        //   gridIndex: 1,
        //   name: yName1,
        //   nameLocation: yInv ? 'start' : 'end',
        //   nameTextStyle: {
        //     fontSize: 16,
        //   },
        //   axisLine: {
        //     onZero: false,
        //   },
        //   axisLabel: {
        //     fontSize: 14,
        //   },
        //   position: 'left',
        //   inverse: yInv,
        //   min: Math.ceil(Math.min.apply(null, embd1_y) - 1.5),
        //   max: Math.ceil(Math.max.apply(null, embd1_y) + 1.5),
        // },
      ],
    }
  }


  useImperativeHandle(onRef, () => ({  // explode trigger for parent components
    "Trigger": toggleAnno, // Trigger for useEffect
    "Loading": setLoadings, // set Loading status
    "Tip": setCurrTip, // set Loading Tips
  }))

  useEffect(() => {
    if (Init) {
      var myChart = echarts.getInstanceByDom(chartRef.current)
      var myChart1 = echarts.getInstanceByDom(chartRef1.current)
      var myChart2 = echarts.getInstanceByDom(chartRef1.current)

      if (commandRef.current === "Reload") {
        enterLoading(0, setLoadings)
        setCurrTip(loadingTips[1])
        setCompLoad((compLoad) => {
          let newCompLoad = { ...compLoad }
          newCompLoad['NetWorkRelation'] = true
          return newCompLoad
        })
        ScH5adLoader(scfile, false).then((res) => {
          console.log(res)
          setSCGraph(res)
          setCellArray(res.categories)
        })
        ScH5adLoader(spfile, true).then((res) => {
          console.log("All data reloaded in Network.")
          setSPGraph(res)
          setGraph(res)
          toggleAnno("Upload")
        })

        // if(key){
        //   myChart1.getDom().style.display = 'block';
        //   myChart2.getDom().style.display = 'none';
        // }
        // else{
        //   myChart1.getDom().style.display = 'none';
        //   myChart2.getDom().style.display = 'block';
        // }
      }
      if (commandRef.current === "Upload") {
        enterLoading(0, setLoadings)

        console.log(graph)

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
                edgeLength: [500, 200],
                layoutAnimation: false
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

        console.log(cellArray)


        // pair view
        if (key) {
          // 1.set source
          let _dims = [...Object.keys(_data[0]), 'id']
          let source = _data.map((item, id) => {
            return [...Object.entries(item).map(([_, value]) => value), id]
          })
          symbolSizeRef.current = source.length > 5000 ? 4 : 6
          setItemSize(symbolSizeRef.current)
          setCellNum(source.length)

          console.log(_dims)
          console.log(graph)

          // 2.0 set annotations and batches
          let defaultAnno = 'annotation' // 'annotation' as annotation
          setClusterCur({ value: 0, label: defaultAnno, attr: 'categories' })
          prevCluster.current = { value: 0, label: defaultAnno, attr: 'categories' }
          // the result saved in itemGroupRef.current
          let annotations = setItemGroup(source, _dims.indexOf(defaultAnno))
          console.log(annotations)
          console.log(itemGroupRef.current)
          let clusterGroups = itemGroupRef.current

          // let batches = ['batch 1']
          // let _batchName = null
          // setBatchName(null)
          // setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
          // setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))

          let batches = []
          let _batchName = batchName
          if (_dims.includes(_batchName)) {
            batches = setItemGroup(source, _dims.indexOf(_batchName))
            setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
            setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))
          }
          else {
            _batchName = null
            setBatchName(null)
            batches = ['batch 1']
            setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
            setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))
          }

          console.log(batchCur)
          console.log(batches)

          // 3. set embeddings
          let _embd = embedOps.map((item) => item.label).includes('spatial')
            ? embedOps.find((item) => item.label === 'spatial')
            : embedOps[0]
          let xName = `${_embd.label}_0`
          let yName = `${_embd.label}_1`
          let axis = setAxis(source, _dims, xName, yName)
          setEmbedCur(_embd)

          // 4. set datasets
          let _datasets = setBatchDataset(source, _dims, defaultAnno, annotations, batches[0], _batchName)

          // 5. set anno series with batch
          let _series = []
          let _snum = 0
          for (let anno in annotations) {
            _snum = _snum + 1
            _series.push({
              type: 'scatter',
              symbolSize: symbolSizeRef.current,
              xAxisIndex: 0,
              yAxisIndex: 0,
              name: anno,
              encode: {
                x: xName,
                y: yName,
                tooltip: [0, 1, 2, 3],
                itemName: anno,
              },
              emphasis: {
                focus: 'series',
              },
              itemStyle: {
                // opacity: itemOpacity,
              },
              large: true,
              largeThreshold: 5000,
              datasetIndex: _snum,
            })
          }

          console.log(axis)
          console.log(_datasets)
          console.log(_series)

          // 6.calculate top3 celltypes
          var celltypeIndexDic = {}
          cellArray.forEach((item) => {
            // find the index in dims
            const index = _dims.indexOf(item.name)
            if(index != -1){
              celltypeIndexDic[item.name] = index
            }
          })
          var cellScores = [] // cellScore[0] = {"CellA":1,"CellB":2.4....}
          // itemGroupRef.current: [[1,3,4,..],[2,6,..]]
          for(let clusterNum = 0;clusterNum < clusterGroups.length;clusterNum++){
            var tempGroup = clusterGroups[clusterNum].map(index => source[index])
            cellScores[clusterNum] = {}
            for(let cellType in celltypeIndexDic){
              let cellIndex = celltypeIndexDic[cellType]
              // const sum = data.reduce((acc, current) => acc + current[1], 0);
              let score = tempGroup.reduce((acc, current) => acc + current[cellIndex], 0) / tempGroup.length
              cellScores[clusterNum][cellType] = score
            }
          }
          console.log(cellScores)
          var topCells = cellScores.map(scores => {
            return Object.entries(scores)
              .sort(([, valueA], [, valueB]) => valueB - valueA)  // sort
              .slice(0, 3)  // top 3
              .map(([key]) => key);  // cell name
          });
          console.log(topCells);



          // 7.add markers
          let _data2 = []
          let showTop = 5
          for (let i = 0; i < cluLen; i++) {
            let _arr = mn.slice(i * geneNum, i * geneNum + showTop)  // Top marker for each cluster
            let _obj = _arr.map((item, index) => ({
              'name': item,
              'lfc': lfc[i * geneNum + index],
              'pval': pa[i * geneNum + index],
              'expr': me[i * geneNum + index],
              'frac': fr[i * geneNum + index],
              'clu': i.toString(),
            }))
            _data2 = [..._data2, ..._obj]
          }
          var topMarkers = []
          console.log(annotations)
          for(let i=0;i<annotations.length;i++){
            topMarkers[i] = []
          }
          for(let i = 0; i < _data2.length; i++){
            topMarkers[Math.floor(i/5)][i%5] = _data2[i].name
          }
          console.log(topMarkers)

          var newToolTip = {
            formatter: function (params) {
              const index = params.dataIndex; 
              const markers = topMarkers[index].join(', ');  // 拼接 topMarkers
              const cells = topCells[index].join(', ');      // 拼接 topCells
              
              const marker = params.marker;
          
              return `${marker} ${index}<br/>` +
                     `<b> Top Markers: </b> ${markers}<br/>` +
                     `<b> Top Cells: </b> ${cells}`;
            }
          }
          let option = myChart.getOption()
          option.tooltip = newToolTip
          myChart.setOption(option)


          // set options
          // let axis = Axis.setEmptyAxis(0)
          myChart1.setOption({
            // title: [
            //   {
            //     text: title,
            //   }
            // ],
            grid: [{
              // top: '14%',
              height: '60%',
              width: '60%',
              // right: '1%',
              left: '10%',
              // bottom: '18%',
            }],
            // backgroundColor: 'green',
            xAxis: axis.xAxis,
            yAxis: axis.yAxis,
            dataset: _datasets,
            series: _series,
          })

          // hover function
          myChart.on('click', (params) => {
            let _seriesName = params.name
            console.log(_seriesName)
            let option = myChart1.getOption()
            option.series.forEach(seriesItem => {
              if (seriesItem.name === _seriesName) {
                seriesItem.itemStyle = {
                  // color: 'red',
                  opacity: 1,
                }
              }
              else {
                seriesItem.itemStyle = {
                  opacity: 0.1
                }
              }

            })
            myChart1.setOption(option)


            let option1 = myChart.getOption()
            option1.series[0].links.forEach(link => {
              if (link.source === _seriesName || link.target === _seriesName) {
                link.lineStyle = {
                  ...link.lineStyle,
                  opacity: 1 // 选中节点相关的边保持正常显示
                }
              } else {
                link.lineStyle = {
                  ...link.lineStyle,
                  opacity: 0.1 // 不相关的边透明度降低
                }
              }
            })
            myChart.setOption(option1)
          })

          // // hover out function
          // myChart.on('mouseout', () => {
          //   let option = myChart1.getOption();
          //   option.series.forEach(seriesItem => {
          //     seriesItem.itemStyle = {
          //       opacity: 1,
          //     };
          //   });
          //   myChart1.setOption(option);
          // })



          // mychart1 hover
          myChart1.on('click', (params) => {
            console.log(graph.nodes)
            console.log(params)
            let option = myChart.getOption()
            let selectedNodes = option.series[0].data.filter(node => node.category === params.seriesName)
            console.log(selectedNodes)
            // selectedNodes.forEach(node => {
            //   node.itemStyle = {

            //   };
            // });
            let selectedNode = selectedNodes[0]

            option.series[0].links.forEach(link => {
              if (link.source === selectedNode.name || link.target === selectedNode.name) {
                link.lineStyle = {
                  ...link.lineStyle,
                  opacity: 1 // 选中节点相关的边保持正常显示
                }
              } else {
                link.lineStyle = {
                  ...link.lineStyle,
                  opacity: 0.1 // 不相关的边透明度降低
                }
              }
            })

            myChart.setOption(option)

            let option1 = myChart1.getOption()
            option1.series.forEach(seriesItem => {
              if (seriesItem.name === selectedNode.name) {
                seriesItem.itemStyle = {
                  // color: 'red',
                  opacity: 1,
                }
              }
              else {
                seriesItem.itemStyle = {
                  opacity: 0.1
                }
              }
            })
            myChart1.setOption(option1)

          })


          // myChart1.on('mouseout', (params) => {
          //   // console.log(graph.nodes)
          //   let option = myChart.getOption();

          //   option.series[0].links.forEach(link => {
          //     link.lineStyle = {
          //       ...link.lineStyle,
          //       opacity: 1 
          //     };
          //   });

          //   myChart.setOption(option);
          // })
        }
        else {
          // lasso view
          // setting for mychart2

          let _dims = [...Object.keys(_data1[0]), 'id']
          let source = _data1.map((item, id) => {
            return [...Object.entries(item).map(([_, value]) => value), id]
          })
          setCellNum(source.length)
          symbolSizeRef1.current = source.length > 5000 ? 2 : 4
          setItemSize1(symbolSizeRef1.current)
          let defaultAnno = clusterOps.find((item) => item.label === "annotation")
          if (typeof defaultAnno === 'undefined') {
            defaultAnno = clusterOps[0]
          }
          setClusterCur(defaultAnno)
          let annotations = setItemGroup(source, _dims.indexOf(defaultAnno.label))

          // set batches
          let batches = []
          if (_dims.includes(batchName)) {
            batches = setItemGroup(source, _dims.indexOf(batchName), 'categories', false)
            setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
            setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))
          }
          else {
            setBatchName(null)
            batches = ['batch 1']
            setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
            setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))
          }

          // set embeddings and axis
          let defaultEmbd = embedOps.map((item) => item.label).includes('X_umap')
            ? 'X_umap'
            : embedOps[0].label
          setEmbedCur(defaultEmbd)
          setEmbedCur1(defaultEmbd)
          let xName = `${defaultEmbd}_0`
          let yName = `${defaultEmbd}_1`
          let axis = setAxis1(source, _dims, xName, yName, xName, yName)

          // set Dataset
          let _datasets = setDataset(
            source,
            _dims,
            defaultAnno.label,
            annotations
          )

          // set Series
          let _series = []
          let _snum = 0
          // _series.push({
          //   type: 'scatter',
          //   xAxisIndex: 1,
          //   yAxisIndex: 1,
          //   encode: {
          //     x: xName,
          //     y: yName,
          //   },
          //   tooltip: { show: false },
          //   itemStyle: {
          //     color: 'gray',
          //     opacity: 0.3,
          //   },
          //   large: true,
          //   largeThreshold: 0,
          //   datasetIndex: _snum,
          // })
          for (let anno of annotations) {
            // seperate the _data into subgroups by cell-types
            _snum = _snum + 1
            _series.push({
              type: 'scatter',
              symbolSize: symbolSizeRef1.current,
              name: anno,
              encode: {
                x: xName,
                y: yName,
                tooltip: [0, 1, 2, 3],
                itemName: anno,
              },
              emphasis: {
                focus: 'series',
              },
              large: true,
              largeThreshold: 5000,
              datasetIndex: _snum,
            })
          }

          // top markers
          let _data2 = []
          let showTop = 10
          for (let i = 0; i < cluLen1; i++) {
            let _arr = mn1.slice(i * geneNum, i * geneNum + showTop)  // Top marker for each cluster
            let _obj = _arr.map((item, index) => ({
              'name': item,
              'lfc': lfc1[i * geneNum + index],
              'pval': pa1[i * geneNum + index],
              'expr': me1[i * geneNum + index],
              'frac': fr1[i * geneNum + index],
              'clu': i.toString(),
            }))
            _data2 = [..._data2, ..._obj]
          }
          console.log(_data2)
          var topMarkers = []
          console.log(annotations)
          for(let i=0;i<annotations.length;i++){
            topMarkers[i] = []
          }
          for(let i = 0; i < _data2.length; i++){
            topMarkers[Math.floor(i/showTop)][i%showTop] = _data2[i].name
          }
          console.log(topMarkers)

          var newToolTip = {
            formatter: function (params) {
              let cellName = params.name
              const index = cellArray.findIndex(cell => cell.name === cellName)
              // const markers = topMarkers[index].join(', ');  // 拼接 topMarkers
              const markers = topMarkers[index];  // 拼接 topMarkers              
              const marker = params.marker;
              const firstHalf = markers.slice(0, 4).join(', ') +',';
              const secondHalf = markers.slice(6).join(', ');
          
              return `${marker} ${cellName}<br/>` +
                     `<b> Top Markers: </b> ${firstHalf}<br/>
                      &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp
                     ${secondHalf}`; 
            }
          }
          let option = myChart.getOption()
          option.tooltip = newToolTip
          myChart.setOption(option)


          myChart2.setOption(
            {
              // title: [
              //   {
              //     text: title,
              //     left: 'center',
              //     textStyle: {
              //       fontSize: 24,
              //     },
              //   },
              // ],
              xAxis: axis.xAxis,
              yAxis: axis.yAxis,
              dataset: _datasets,
              series: _series,
              // toolbox: {
              //   top: yInv ? '4%' : '7%',
              // },
            }
          )

          // hover function
          myChart.on('click', (params) => {
            let _seriesName = params.name
            let option = myChart2.getOption()
            option.series.forEach(seriesItem => {
              if (seriesItem.name === _seriesName) {
                seriesItem.itemStyle = {
                  // color: 'red',
                  opacity: 1,
                }
              }
              else {
                seriesItem.itemStyle = {
                  opacity: 0.1
                }
              }

            })
            myChart2.setOption(option)
          })

          // // hover out function
          // myChart.on('mouseout', () => {
          //   let option = myChart2.getOption();
          //   option.series.forEach(seriesItem => {
          //     seriesItem.itemStyle = {
          //       opacity: 1,
          //     };
          //   });
          //   myChart2.setOption(option);
          // })


          // mychart2 hover
          myChart2.on('click', (params) => {
            console.log("mychart2", params)
            // console.log(graph.nodes)
            let option = myChart.getOption()
            let selectedNodes = option.series[0].data.filter(node => node.category === params.seriesName)
            console.log(selectedNodes)
            // selectedNodes.forEach(node => {
            //   node.itemStyle = {

            //   };
            // });
            let selectedNode = selectedNodes[0]

            option.series[0].links.forEach(link => {
              if (link.source === selectedNode.name || link.target === selectedNode.name) {
                link.lineStyle = {
                  ...link.lineStyle,
                  opacity: 1 // 选中节点相关的边保持正常显示
                }
              } else {
                link.lineStyle = {
                  ...link.lineStyle,
                  opacity: 0.1 // 不相关的边透明度降低
                }
              }
            })

            myChart.setOption(option)
          })

          // myChart2.on('mouseout', (params) => {
          //   let option = myChart.getOption();

          //   option.series[0].links.forEach(link => {
          //     link.lineStyle = {
          //       ...link.lineStyle,
          //       opacity: 1 
          //     };
          //   });

          //   myChart.setOption(option);
          // })

        }

        // if(key){
        //   myChart1.getDom().style.display = 'block';
        //   myChart2.getDom().style.display = 'none';
        // }
        // else{
        //   myChart1.getDom().style.display = 'none';
        //   myChart2.getDom().style.display = 'block';
        // }

        quitLoading(0, setLoadings)
      }
      else if (commandRef.current == "Setting") {
        if (key) {
          let option = myChart1.getOption()
          console.log(option)
          let _source = option.dataset[0].source
          let _dims = option.dataset[0].dimensions
          let bat = batchCur.label
          let clu = clusterCur.label
          let cluIdx = _dims.indexOf(clu)
          let embd = embedCur.label
          let _snum = 1
          let _series = option.series
          let xName = `${embd}_0`
          let yName = `${embd}_1`
          let axis = setAxis(_source, _dims, xName, yName)
          // set Annotations
          let annotations = setItemGroup(_source, _dims.indexOf(clu), clusterCur.attr)

          // set Dataset
          let _dataset = setBatchDataset(_source, _dims, clu, annotations, bat, batchName)
          let _newSeries = []
          // setting anno with batch in 2D
          for (let anno in annotations) {
            _newSeries.push({
              type: 'scatter',
              symbolSize: itemSize,
              xAxisIndex: 0,
              yAxisIndex: 0,
              name: anno,
              encode: {
                x: xName,
                y: yName,
                itemName: anno,
              },
              emphasis: {
                focus: 'series',
              },
              itemStyle: {
              },
              large: true,
              largeThreshold: 5000,
              datasetIndex: _snum,
            })
            _snum = _snum + 1
          }

          _series = _newSeries

          myChart1.setOption({
            // title: [
            //   {
            //     text: title,
            //     top: '0%',
            //     left: 'center',
            //   },
            // ],
            xAxis: axis.xAxis,
            yAxis: axis.yAxis,
            series: _series,
            dataset: _dataset,
          },)
        }
        else {
          // lasso settings

          let option = myChart2.getOption()
          let _source = option.dataset[0].source
          let _dims = option.dataset[0].dimensions
          let clusterIdx = _dims.indexOf(clusterCur.label)

          // set axis
          let xName0 = embedCur + '_0'
          let yName0 = embedCur + '_1'
          let xName1 = embedCur1 + '_0'
          let yName1 = embedCur1 + '_1'
          let axis = setAxis(_source, _dims, xName0, yName0, xName1, yName1)


          let _len = option.series.length
          let defAnnolen = itemGroupRef.current.length
          let _series = Array(_len)
            .fill()
            .map((_, id) => ({
              symbolSize: itemSize1,
              // itemStyle: { ...itemStyle},
              encode: {
                x: id > 0 && id <= defAnnolen ? xName0 : xName1,
                y: id > 0 && id <= defAnnolen ? yName0 : yName1,
              },
            }))
          myChart2.setOption({
            xAxis: axis.xAxis,
            yAxis: axis.yAxis,
            series: _series,
          })

        }



      }
    } else {
      var myChart = echarts.init(chartRef.current) //init the echart container
      var myChart1 = echarts.init(chartRef1.current)
      var myChart2 = echarts.init(chartRef1.current)

      enterLoading(0, setLoadings)
      setCurrTip(loadingTips[0])
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
        // title: [
        //   {
        //     text: title,
        //     left: '90%',
        //     textStyle: {
        //       fontSize: 24,
        //     },
        //   },
        // ],
        animation: false,
        grid: [{
          top: '0%',
          bottom: '0%',
        }],
        color: vega_20,
      })

      let axis = Axis.setEmptyAxis(0)
      myChart1.setOption({
        textStyle: {
          fontFamily: 'Arial'
        },
        tooltip: {},
        animation: false,
        // title: [
        //   {
        //     text: title,
        //     top: '0%',
        //     left: 'center',
        //   },
        // ],
        xAxis: axis.xAxis,
        yAxis: axis.yAxis,
        grid: [{
          top: '18%',
          width: '80%',
          right: '1%',
          bottom: '18%',
        }],
        color: vega_20

      })

      let axis0 = Axis.setEmptyAxis(0)
      let axis1 = Axis.setEmptyAxis(1)

      myChart2.setOption({
        tooltip: {},
        // brush: {
        //   toolbox: ['rect', 'polygon', 'clear'],
        //   outOfBrush: {
        //     colorAlpha: 0.3,
        //   },
        //   xAxisIndex: 0,
        //   throttleType: 'debounce',
        //   throttleDelay: 1000,
        // },
        // toolbox: {
        //   show: true,
        //   itemSize: 20,
        //   itemGap: 10,
        //   feature: {
        //     bursh: {
        //       type: ['rect', 'polygon', 'clear'],
        //     },
        //     mark: { show: true },
        //     dataView: { show: true, readOnly: true },
        //     restore: { show: true },
        //     saveAsImage: { show: true, pixelRatio: 5, },
        //     dataZoom: {},
        //   },
        //   iconStyle: {
        //     borderWidth: 1.5,
        //   },
        //   top: '7%',
        // },
        // title: [
        //   {
        //     text: title,
        //     left: 'center',
        //     textStyle: {
        //       fontSize: 24,
        //     },
        //   },
        // ],
        grid: [
          {
            top: '18%',
            left: '5%',
            width: '43%',
            bottom: '18%',
          },
          {
            top: '18%',
            width: '43%',
            right: '1%',
            bottom: '18%',
          },
        ],
        // legend: {
        //   type: 'scroll',
        //   animation: false,
        //   pageIconColor: '#4096ff',
        //   orient: 'horizontal',
        //   bottom: '1%',
        //   padding: 0,
        //   textStyle: {
        //     fontSize: 16,
        //   },
        // },
        color: vega_20
      })



      setInit(true)
    }
  }, [action, key])




  return (
    <Spin
      spinning={loadings[0]}
      size="large"
      tip={
        <div>{currTip}
          <Progress percent={progress} strokeColor={strokeColor} size={[300, 15]} />
        </div>
      }>
      <div>
          <h2 style={{ fontSize: '24px', marginLeft: '15rem' }}>{"Cell-Cell Interaction Network"}</h2>
      </div>
      <div style={{ display: 'flex', width: '80rem'}}>
        <div
          ref={chartRef}
          className="chart"
          //the target DOM container needs height and width
          style={{ height: height, width: '40rem', marginLeft: '-15rem' }}>
        </div>
        <div
          ref={chartRef1}
          className="chart1"
          style={{ height: height, width: '40rem' }}>
        </div>
        <div style={{ 
        // display: 'flex' 
          // marginRight: '10rem',
          // width: '20rem',
          display: 'flex', flexDirection: 'column', marginLeft: '-8rem',
        }}>
        <Space direction='vertical'>
          <Switch checkedChildren="Spatial"
            unCheckedChildren="Single-cell"
            style={{
              // marginRight:'5rem',
              marginTop: '20rem'
            }}
            value={key}
            onChange={(value) => {
              // recovery
              var myChart1 = echarts.getInstanceByDom(chartRef1.current)
              let option = myChart1.getOption()
              option.series.forEach(seriesItem => {
                seriesItem.itemStyle = {
                  opacity: 1,
                }
              })
              myChart1.setOption(option)

              var myChart = echarts.getInstanceByDom(chartRef.current)
              let myoption = myChart.getOption()
              myoption.series[0].links.forEach(link => {
                link.lineStyle = {
                  ...link.lineStyle,
                  opacity: 1 // 选中节点相关的边保持正常显示
                }
              })
              myChart.setOption(myoption)

              // exchange
              setKey(value)
              setGraph(value ? spgraph : scgraph)
              toggleAnno("Upload")
            }} />
          <Button
            type="primary"
            style={{ marginTop: '1rem' }}
            onClick={() => {
              var myChart1 = echarts.getInstanceByDom(chartRef1.current)
              let option = myChart1.getOption()
              option.series.forEach(seriesItem => {
                seriesItem.itemStyle = {
                  opacity: 1,
                }
              })
              myChart1.setOption(option)

              var myChart = echarts.getInstanceByDom(chartRef.current)
              let myoption = myChart.getOption()
              myoption.series[0].links.forEach(link => {
                link.lineStyle = {
                  ...link.lineStyle,
                  opacity: 1 // 选中节点相关的边保持正常显示
                }
              })
              myChart.setOption(myoption)
            }}
          >
            {"recovery"}
          </Button>
          <div style={{ marginTop: '0.5rem'}}
          >Batch:</div>
          <Select
            type="primary"
            labelInValue
            value={batchCur.label}
            options={batchOps}
            onChange={(target) => {
              // console.log(target)
              let selectedIndex = target.value
              setBatchCur(batchOps[selectedIndex])
              toggleAnno("Setting")
            }}
          ></Select>
      </Space>
      </div>
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
  setCompLoad: PropTypes.func,
  onRef: PropTypes.any,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
  progress: PropTypes.number,
}

export default NetworkRelation

