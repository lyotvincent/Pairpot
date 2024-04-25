import React, { useRef, useEffect, useState } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import { TooltipComponent, VisualMapComponent } from 'echarts/components'
import {
  ConfigProvider,
  Button,
  Card,
  Space,
  Popconfirm,
  Input,
  Select,
  notification,
  Popover,
  Upload,
  Segmented,
  InputNumber,
  Slider,
  Switch,
  Col,
  Row,
  Spin,
  Form,
  Flex
} from 'antd'
import {
  GraphicComponent,
  GridComponent,
  LegendComponent,
} from 'echarts/components'
import { LineChart } from 'echarts/charts'
import { UniversalTransition } from 'echarts/features'
import { CanvasRenderer } from 'echarts/renderers'
import H5adLoader from '../utils/H5adLoader'
import { CloudDownloadOutlined, FileAddOutlined, FileOutlined, InboxOutlined, SettingOutlined, ReloadOutlined } from '@ant-design/icons'
echarts.use([
  GraphicComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  LineChart,
  CanvasRenderer,
  UniversalTransition,
])

const { Dragger } = Upload

// left figure is umap2 and right figure is spatial3
const SpScatter = ({ theme, height, width, margin }) => {
  const chartRef = useRef(null) // get current DOM container
  const commandRef = useRef('') // get the current command for useEffect
  const [action, setAction] = useState(0)
  const [api, contextHolder] = notification.useNotification()
  const [isInit, setInit] = useState(false) // whether echart object is inited
  const [clusterCur, setClusterCur] = useState({})
  const [batchName, setBatchName] = useState("batch")
  const [batchCur, setBatchCur] = useState({})
  const [embedCur, setEmbedCur] = useState('')
  const [embed3Dcur, setEmbed3Dcur] = useState('')
  const [clusterOps, setClusterOps] = useState([])
  const [embedOps, setEmbedOps] = useState([])
  const [batchOps, setBatchOps] = useState([])
  const [loading, setLoading] = useState(false)
  const [itemSize, setItemSize] = useState(2)
  const [itemOpacity, setItemOpacity] = useState(0.8)
  const [xInv, setxInv] = useState(false)
  const [yInv, setyInv] = useState(false)
  const [cellNum, setCellNum] = useState(0)
  const itemGroupRef = useRef([])
  const symbolSizeRef = useRef('')
  const FileLoaderRef = useRef('')
  const prevCluster = useRef('')
  const [title, setTitle] = useState('Mouse-Brain')
  const [_data, _setData] = useState(
    require('../../assets/data/Mouse-Brain-umap.json')
  )
  var vega_20 = [
    '#1f77b4',
    '#aec7e8',
    '#ff7f0e',
    '#ffbb78',
    '#2ca02c',
    '#98df8a',
    '#d62728',
    '#ff9896',
    '#9467bd',
    '#c5b0d5',
    '#8c564b',
    '#c49c94',
    '#e377c2',
    '#f7b6d2',
    '#7f7f7f',
    '#c7c7c7',
    '#bcbd22',
    '#dbdb8d',
    '#17becf',
    '#9edae5',
  ]

  vega_20 = [...vega_20, ...vega_20]
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
      xAxis3D: [{
        grid3DIndex: 0,
        name: xName0,
        nameLocation: 'middle',
        min: Math.ceil(Math.min.apply(null, embd0_x) - 1.5),
        max: Math.ceil(Math.max.apply(null, embd0_x) + 1.5),
      }
      ],
      yAxis3D: [{
        grid3DIndex: 0,
        name: yName0,
        nameLocation: 'middle',
        min: Math.ceil(Math.min.apply(null, embd0_y) - 1.5),
        max: Math.ceil(Math.max.apply(null, embd0_y) + 1.5),
      }
      ],
    }
  }

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

  const setDataset = (source, dims, annoName, annotations) => {
    let _datasets = []
    _datasets.push({
      dimensions: dims,
      source: source,
    })
    if (annotations.length === 0) {
      // if no annotations, all the data are annotated by a label.
      _datasets.push({
        // 这个 dataset 的 index 是 `1`。
        transform: {
          type: 'sort',
          config: { dimension: annoName, order: 'desc' },
        },
      })
    } else {
      for (let anno of annotations) {
        _datasets.push({
          // 这个 dataset 的 index 是 `1`。
          transform: {
            type: 'filter',
            config: { dimension: annoName, value: anno },
          },
        })
      }
    }
    return _datasets
  }

  const setBatchDataset = (source, dims, annoName, annotations, batch, batchName) => {
    let _datasets = []
    _datasets.push({
      dimensions: dims,
      source: source,
    })
    if(batchName !== null){
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
    }else{
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

  const SpH5adLoader = (file) => {
    const reader = new FileReader()
    reader.onload = (event) => {
      let h5info = H5adLoader(file, event)
      _setData(h5info.data)
      setTitle(h5info.title)
      setClusterOps(h5info.clusters)
      setEmbedOps(h5info.embdOps)
    }
    reader.onloadend = () => {
      toggleAnno('Upload')
      setLoading(false)
    }
    reader.readAsArrayBuffer(file)
  }

  const beforeUpload = (file) => {
    if (file.name.endsWith('.h5ad')) {
      FileLoaderRef.current = SpH5adLoader
    } else {
      alert('Please select a .json or .h5ad file')
      return false
    }
    return true
  }

  const onUpload = (info) => {
    const file = info.file
    setLoading(true)
    FileLoaderRef.current(file)
  }

  const upLoadProps = {
    maxCount: 1,
    beforeUpload: beforeUpload,
    customRequest: onUpload,
    listType: 'text',
    showUploadList: false,
  }

  const toggleAnno = (command) => {
    commandRef.current = command
    setAction(action + 1)
  }

  useEffect(() => {
    if (isInit) {
      console.log("Is Inited.")
      var myChart = echarts.getInstanceByDom(chartRef.current)
      if (commandRef.current === "Upload") {

        // 1.set source
        let _dims = [...Object.keys(_data[0]), 'id']
        let source = _data.map((item, id) => {
          return [...Object.entries(item).map(([_, value]) => value), id]
        })
        symbolSizeRef.current = source.length > 5000 ? 2 : 4
        setItemSize(symbolSizeRef.current)
        setCellNum(source.length)

        // 2.set annotations and batches
        let defaultAnno = 'annotation'
        setClusterCur({ value: 0, label: defaultAnno, attr: 'categories' })
        prevCluster.current = { value: 0, label: defaultAnno, attr: 'categories' }
        let annotations = setItemGroup(source, _dims.indexOf(defaultAnno))

        let batches = []
        if(_dims.includes(batchName)){
          batches = setItemGroup(source, _dims.indexOf(batchName))
          setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
          setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))
        }
        else {
          setBatchName(null)
          batches = ['batch 1']
          setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
          setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))
        }


        // 3.set embeddings
        let _embd = embedOps.map((item) => item.label).includes('spatial')
        ? 'spatial'
        : embedOps[0]
        let xName = `${_embd.label}_0`
        let yName = `${_embd.label}_1`
        let axis = setAxis(source, _dims, xName, yName)
        setEmbedCur(_embd)
        setEmbed3Dcur(_embd)

        // 4.set datasets
        let _datasets = setBatchDataset(source, _dims, defaultAnno, annotations, batches[0], batchName)
        console.log(_datasets)
        // 5.set 3D series
        let _series = []
        let _snum = 0
        let _3DSeries = {
          type: 'scatter3D',
          coordinateSystem: 'cartesian3D',
          symbolSize: symbolSizeRef.current,
          name: "batch",
          encode: {
            x: xName,
            y: yName,
            z: batchName,
            tooltip: [0],
          },
          itemStyle: {
            color: 'gray',
            opacity: 0.8,
          },
          datasetIndex: _snum,
          zlevel: -1,
        }
        _series.push(_3DSeries)
        // 6.set anno series with batch
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
              opacity: itemOpacity,
            },
            large: true,
            largeThreshold: 5000,
            datasetIndex: _snum,
          })
        }

        // 7.set anno series without batch
        for (let anno in annotations) {
          _snum = _snum + 1
          _series.push({
            type: 'scatter3D',
            symbolSize: symbolSizeRef.current,
            name: anno,
            encode: {
              x: xName,
              y: yName,
              z: 'batch',
              tooltip: [0],
              itemName: anno,
            },
            itemStyle: {
              opacity: 0.7,
            },
            emphasis: {
              focus: 'series',
            },
            datasetIndex: _snum,
          })
        }

        // //8. set visualMap
        // let _visualMap = {
        //   type: 'piecewise',
        //   categories: batches,
        //   dimension: _dims.indexOf(defaultBatch),
        //   orient: 'horizontal',
        //   top: "5%",
        //   left: "0%",
        //   right: '50%',
        //   selectedMode: "single",
        //   inRange: {
        //     color: vega_20,
        //   },
        //   outOfRange: {
        //     color: '#808080'
        //   },
        //   seriesIndex: [...Array(annotations).keys()],
        // }


        // 9.render charts
        myChart.setOption({
          title: [
            {
              text: title,
              top: '0%',
              left: 'center',
              textStyle: {
                fontSize: 24,
              },
            },
          ],
          grid3D: {
            top: 'top',
            width: '55%',
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
          xAxis3D: axis.xAxis3D,
          yAxis3D: axis.yAxis3D,
          zAxis3D: { type: 'category' },
          xAxis: axis.xAxis,
          yAxis: axis.yAxis,
          dataset: _datasets,
          series: _series,
        },
          {
            replaceMerge: ['dataset', 'series', 'visualMap'],  // enable replaceMerge for datasets
          })

        // register click action
        myChart.on('click', (params) => {
          if (params.componentType === 'series') {
            let _seriesName = params._seriesName
            myChart.setOption({
              series: {
                name: _seriesName,
                itemStyle: {
                  color: 'red',
                  opacity: 1
                }
              }
            })
          }
        })
      }
      else if (commandRef.current === "Setting") {
        let option = myChart.getOption()
        let _source = option.dataset[0].source
        let _dims = option.dataset[0].dimensions
        let bat = batchCur.label
        let clu = clusterCur.label
        let cluIdx = _dims.indexOf(clu)
        let embd = embedCur.label
        let embd3D = embed3Dcur.label
        let _snum = 1
        let _series = option.series
        let xName = `${embd}_0`
        let yName = `${embd}_1`
        let axis = setAxis(_source, _dims, xName, yName)

        // set Annotations
        let annotations = setItemGroup(_source, _dims.indexOf(clu), clusterCur.attr)

        // set Dataset
        let _dataset = setBatchDataset(_source, _dims, clu, annotations, bat, batchName)
        console.log(_dataset)
        // set series
        let myVisualMap = []
        let _3DSeries = _series[0]
        _3DSeries.symbolSize = itemSize
        _3DSeries.itemStyle.opacity = itemOpacity
        _3DSeries.encode= {
          x: `${embd3D}_0`,
          y: `${embd3D}_1`,
          z: batchName,
          tooltip: [0],
        }
        let _newSeries = []
        if (clusterCur.attr === 'categories') {
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
                opacity: itemOpacity
              },
              large: true,
              largeThreshold: 5000,
              datasetIndex: _snum,
            })
            _snum = _snum + 1
          }
          // setting anno without batch in 3D
          for (let anno in annotations) {
            _newSeries.push({
              type: 'scatter3D',
              symbolSize: itemSize,
              grid3DIndex: 0,
              name: anno,
              encode: {
                x: `${embd3D}_0`,
                y: `${embd3D}_1`,
                z: batchName,
                tooltip: [0],
                itemName: anno,
              },
              itemStyle: {
                opacity: itemOpacity
              },
              datasetIndex: _snum,
            })
            _snum = _snum + 1
          }
        }
        else {
          // add 2D series
          _newSeries.push({
            type: 'scatter',
            name: clusterCur.label,
            symbolSize: itemSize,
            encode: {
              // annotations are displayed in left chart
              x: xName,
              y: yName,
              tooltip: [
                0,
                _dims.indexOf(xName),
                _dims.indexOf(yName),
              ],
            },
            itemStyle: {
              opacity: itemOpacity,
            },
            datasetIndex: _snum,
          })

          // add 3D series
          _newSeries.push({
            type: 'scatter3D',
            name: clusterCur.label,
            symbolSize: itemSize,
            grid3DIndex: 0,
            encode: {
              // annotations are displayed in left chart
              x: `${embd3D}_0`,
              y: `${embd3D}_1`,
              z: batchName,
              tooltip: [0],
            },
            itemStyle: {
              opacity: itemOpacity,
            },
            datasetIndex: 0,
          })
          myVisualMap.push({
            show: true,
            id: clusterCur.label + '_visual',
            calculable: true,
            dimension: clusterCur.label,
            seriesIndex: [1, 2],
            left: '0%',
            top: '10%',
            orient: 'horizontal',
            precision: 2,
            calculable: true,
            min: Math.min.apply(
              null,
              _source.map((item) => item[cluIdx])
            ),
            max: Math.max.apply(
              null,
              _source.map((item) => item[cluIdx])
            ),
            inRange: {
              color: ['#808080', '#FF3300'],
              opacity: [0.8, 1],
            },
            text: [clusterCur.label, ''],
            textGap: 20,
            textStyle: {
              fontSize: 16,
            },
          })
        }
        _series = [_3DSeries, ..._newSeries]

        myChart.setOption({
          title: [
            {
              text: title,
              top: '0%',
              left: 'center',
              textStyle: {
                fontSize: 24,
              },
            },
          ],
          xAxis: axis.xAxis,
          yAxis: axis.yAxis,
          series: _series,
          dataset: _dataset,
          visualMap: myVisualMap,
        },
          {
            replaceMerge: ['series', 'dataset', 'visualMap'],
          })
        console.log(myChart.getOption())
      }
    } else {
      var myChart = echarts.init(chartRef.current, theme) //init the echart container

      let _dims = [...Object.keys(_data[0]), 'id']
      let source = _data.map((item, id) => {
        return [...Object.entries(item).map(([_, value]) => value), id]
      })
      symbolSizeRef.current = source.length > 5000 ? 2 : 4
      setCellNum(source.length)
      // set annotations
      let defaultAnno = 'leiden'
      setClusterCur({ value: 0, label: defaultAnno, attr: 'categories' })
      prevCluster.current = { value: 0, label: defaultAnno, attr: 'categories' }
      setClusterOps([{ value: 0, label: defaultAnno, attr: 'categories' }])
      let annotations = setItemGroup(source, _dims.indexOf('leiden'))

      let batches = []
      if (_dims.includes(batchName)){
        batches = setItemGroup(source, _dims.indexOf(batchName))
        setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
        setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))
      }else {
        api.warning({
          message: `default batch is not find in .obs or .json file`,
          placement: 'topRight',
        })
        setBatchName(null)
        batches = ['batch 1']
        setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
        setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))
      }



      // set embeddings

      let axis = setAxis(source, _dims, 'array_row', 'array_col')
      setEmbedCur('spatial')
      setEmbedOps([{ value: 0, label: 'spatial' }])

      // set datasets
      let _datasets = setDataset(source, _dims, 'annotation', annotations)

      // set Series
      let _series = []
      let _snum = 0

      let option = {
        tooltip: {},
        grid3D: {
          width: '50%',
        },
        xAxis3D: {},
        yAxis3D: {},
        zAxis3D: {},
        grid: [
          {
            top: '15%',
            left: '5%',
            width: '43%',
            bottom: '13%',
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

      }
      _series.push({
        type: 'scatter',
        xAxisIndex: 1,
        yAxisIndex: 1,
        symbolSize: symbolSizeRef.current,
        encode: {
          x: 'umap_0',
          y: 'umap_1',
          z: 'batch',
        },
        tooltip: { show: false },
        itemStyle: {
          color: 'gray',
          opacity: 0.3,
        },
        large: true,
        largeThreshold: 0,
        datasetIndex: _snum,
      })


      var batchSeries = [
        {
          type: 'scatter3D',
          name: "batch",
          encode: {
            x: 'array_row',
            y: 'array_col',
            z: 'in_tissue',
            tooltip: [0, 1, 2],
          },
          itemStyle: {
            color: 'gray',
            opacity: 0.3,
          },
          datasetIndex: _snum,
        }
      ]

      var umapSeries = [
        {
          type: 'scatter',
          symbolSize: 6,
          name: 'SOView',
          datasetIndex: _snum,
          encode: {
            x: 'array_col',
            y: 'array_row',
            tooltip: [0, 1, 2, 3, 4],
          },
          itemStyle: {
            color: (params) => {
              return (
                'rgba(' +
                params.data[1] * 255 +
                ',' +
                params.data[2] * 255 +
                ',' +
                params.data[3] * 255 +
                ',1)'
              )
            },
          },
        },
      ]
      var legendData = []
      var cluSeries = []
      var item = _data[0]
      for (let key of Object.keys(item)) {
        if (
          [
            'BayesSpace',
            'SEDR',
            'stlearn',
            'Giotto',
            'Seurat',
            'leiden',
            'SpaGCN',
            'annotation',
          ].includes(key)
        ) {
          legendData.push(key)
          cluSeries.push({
            type: 'scatter',
            symbolSize: 6,
            name: key,
            dimensions: ['array_col', 'array_row', 'Cluster'],
            data: _data.map((item) => {
              return [item['array_col'], item['array_row'], item[key]]
            }),
            encode: {
              x: 'array_col',
              y: 'array_row',
              tooltip: [0, 1, 2],
            },
            itemStyle: {
              color: (params) => {
                return vega_20[params.data[2]]
              },
            },
          })
        }
      }
      _series = [...batchSeries, ...cluSeries]
      myChart.setOption({
        tooltip: {},
        title: [
          {
            text: title,
            top: '0%',
            left: 'center',
            textStyle: {
              fontSize: 24,
            },
          },
        ],
        grid3D: {
          top: 'top',
          width: '55%',
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
        xAxis3D: axis.xAxis3D,
        yAxis3D: axis.yAxis3D,
        zAxis3D: { type: 'category' },
        xAxis: axis.xAxis,
        yAxis: axis.yAxis,
        grid: [{
          top: '15%',
          width: '40%',
          right: '1%',
          bottom: '15%',
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
        }],
        legend: {
          type: 'scroll',
          animation: false,
          pageIconColor: '#4096ff',
          orient: 'horizontal',
          bottom: '1%',
          padding: 0,
          textStyle: {
            fontSize: 16,
          },
          selector: [{ type: 'all', title: "All" }, { type: 'inverse', title: 'Inverse' }],
        },
        visualMap: {
          type: 'piecewise',
          categories: annotations,
          dimension: _dims.indexOf(defaultAnno),
          orient: 'horizontal',
          top: "5%",
          left: "0%",
          right: '50%',
          selectedMode: "single",
          inRange: {
            color: vega_20,
          },
          outOfRange: {
            color: '#808080'
          },
          seriesIndex: 0
        },
        toolbox: {
          show: true,
          itemSize: 20,
          itemGap: 10,
          feature: {
            mark: { show: true },
            dataView: { show: true, readOnly: true },
            restore: { show: true },
            saveAsImage: { show: true },
            dataZoom: {},
          },
          iconStyle: {
            borderWidth: 1.5,
          },
          top: '5%',
        },
        dataset: _datasets,
        series: _series,
      })
      setInit(true)
    }
  }, [theme, action])

  return (
    <Flex justify="center" gap='middle'>
      <Spin spinning={loading} size="large">
        <div
          ref={chartRef}
          className="chart"
          //the target DOM container needs height and width
          style={{ height: height, width: width, margin: margin }}></div>
      </Spin>
      <ConfigProvider
        theme={{
          components: {
            Form: {
              itemMarginBottom: '4px',
            },
            Button: {
              borderColorDisabled: 'gray',
              colorTextDisabled: 'gray',
              colorBgContainerDisabled: 'rgba(0.5, 0.5, 0.5, 0.1)',
            },
            Upload: {
              padding: 0,
            },
          },
        }}>
        <Space direction='vertical' size='small'>
          <Space direction='vertical' size='small'>
            <Dragger {...upLoadProps}>
              <p
                className="ant-upload-drag-icon"
                style={{ fontSize: 16, fontFamily: 'Arial', color: 'white', margin: 5 }}>
                <InboxOutlined />
                <br />
                Upload
              </p>
            </Dragger>
            <Popover
              title={<>Setting Figure Options:</>}
              placement="topLeft"
              content={
                <Space size="small" direction="vertical">
                  <Form
                    name="Settings"
                    size="large"
                    labelCol={{
                      span: 6,
                    }}
                    wrapperCol={{
                      span: 24,
                    }}>
                    <Form.Item label="Size">
                      {' '}
                      <Slider
                        min={0}
                        max={10}
                        marks={{ 0: '0', 5: '5', 10: '10' }}
                        onChange={(value) => {
                          setItemSize(value)
                        }}
                        value={typeof itemSize === 'number' ? itemSize : 2}
                      />
                    </Form.Item>
                    <Form.Item label="Opacity">
                      {' '}
                      <Slider
                        min={0}
                        max={1}
                        marks={{ 0: '0', 0.5: '0.5', 1: '1.0' }}
                        step={0.01}
                        onChange={(value) => {
                          setItemOpacity(value)
                        }}
                        value={
                          typeof itemOpacity === 'number' ? itemOpacity : 0.8
                        }
                      />
                    </Form.Item>
                    <Form.Item
                      label="Inverse"
                      labelCol={{
                        span: 6,
                      }}>
                      <Space>
                        xAxis:
                        <Switch
                          checked={xInv}
                          onChange={() => {
                            setxInv(!xInv)
                          }}
                        />
                        yAxis:
                        <Switch
                          checked={yInv}
                          onChange={() => {
                            setyInv(!yInv)
                          }}
                        />
                      </Space>
                    </Form.Item>
                    <Form.Item
                      label="Batch"
                      labelCol={{
                        span: 6,
                      }}>
                      {' '}
                      <Select
                        labelInValue
                        placeholder="Batch"
                        placement="topLeft"
                        options={batchOps}
                        value={batchCur}
                        onChange={(target) => {
                          setBatchCur({
                            value: target.value,
                            label: target.label,
                            attr: batchOps[target.value].attr,
                          })
                          toggleAnno("Setting")
                        }}
                      />
                    </Form.Item>
                    <Form.Item
                      label="Cluster"
                      labelCol={{
                        span: 6,
                      }}>
                      {' '}
                      <Select
                        labelInValue
                        placeholder="Cluster"
                        style={{
                          width: '100%',
                        }}
                        placement="topLeft"
                        options={clusterOps}
                        value={clusterCur}
                        onChange={(target) => {
                          setClusterCur({
                            value: target.value,
                            label: target.label,
                            attr: clusterOps[target.value].attr,
                          })
                          toggleAnno("Setting")
                        }}
                      />
                    </Form.Item>
                    <Form.Item
                      label="Coord"
                      labelCol={{
                        span: 6,
                      }}>
                      {' '}
                      <Select
                        labelInValue
                        placeholder="Coord"
                        style={{
                          width: '100%',
                        }}
                        placement="topLeft"
                        options={embedOps}
                        value={embedCur}
                        onChange={(target) => {
                          setEmbedCur({
                            value: target.value,
                            label: target.label,
                            attr: embedOps[target.value].attr,
                          })
                          toggleAnno("Setting")
                        }}
                      />
                    </Form.Item>
                    <Form.Item
                      label="Coord3D"
                      labelCol={{
                        span: 6,
                      }}>
                      {' '}
                      <Select
                        labelInValue
                        placeholder="Coord3D"
                        style={{
                          width: '100%',
                        }}
                        placement="topLeft"
                        options={embedOps}
                        value={embed3Dcur}
                        onChange={(target) => {
                          setEmbed3Dcur({
                            value: target.value,
                            label: target.label,
                            attr: embedOps[target.value].attr,
                          })
                          toggleAnno("Setting")
                        }}
                      />
                    </Form.Item>
                  </Form>
                  <Space size="small">
                    <Button icon={<ReloadOutlined />}>Reset</Button>
                    <Button
                      type="primary"
                      icon={<SettingOutlined />}
                      onClick={() => {
                        toggleAnno('Setting')
                      }}>
                      Apply
                    </Button>
                  </Space>
                </Space>
              }
              trigger="click">
              <Button type="primary" block icon={<SettingOutlined />}>
                Settings
              </Button>
            </Popover>
            <Button type="primary" block icon={<CloudDownloadOutlined />}>
              Save
            </Button>
          </Space>
          <div>
            <b>Current Status</b>
            <div>Species: {"Mouse"}</div>
            <div>Organs: {"Brain"}</div>
            <div>Batch: {batchCur.label}</div>
            <div>Clustering: {clusterCur.label}</div>
            <div>Embedding: {embedCur.label}</div>
            <div>Attributes: {clusterOps.length}</div>
            <div>Batches: {batchOps.length}</div>
            <div>Spot number: {cellNum}</div>
            <div>Spot Radius: {"55um"}</div>
          </div>
        </Space>
      </ConfigProvider>
    </Flex>
  )
}

SpScatter.defaultProps = {
  theme: 'dark',
  height: '31rem',
  width: '55rem',
  margin: '2rem',
}

SpScatter.propTypes = {
  theme: PropTypes.string,
  title: PropTypes.string,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default SpScatter
