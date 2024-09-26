import React, { useRef, useEffect, useState, useImperativeHandle } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import strokeColor from '../theme/strokeColor'
import { TooltipComponent, VisualMapComponent } from 'echarts/components'
import Axis from './Axis'
import loadingTips from './LoadingTip'
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
  theme,
  Slider,
  Switch,
  Col,
  Row,
  Spin,
  Form,
  Flex,
  Progress,
  Tour
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
import Loading from './Loading'
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
const { useToken } = theme
const { quitLoading, enterLoading } = Loading

// left figure is umap2 and right figure is spatial3
const SpScatter = ({ spfile, setCompLoad, onRef, height, width, margin, meta, progress }) => {
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
  const [loading, setLoading] = useState([])
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
  const [_data, _setData] = useState([])
  const { token } = useToken()
  const [tourOpen, setTourOpen] = useState(false)

  // loading tips
  const [currTip, setCurrTip] = useState(loadingTips[0])

  // ref for tours
  const UploadRef = useRef(null)
  const SettingsRef = useRef(null)
  const SaveRef = useRef(null)
  const StatusRef = useRef(null)
  const steps = [
    {
      title: 'Upload File',
      description: "Click to upload your file, or drag your file here. Currently, this option accepts .h5ad file with basical attributes: .obs[\'annotation\', \'batch\'], and .obsm['X_umap', 'spatial']. These attributes can be obtained by `Scanpy` pipeline.",
      target: () => UploadRef.current,
    },
    {
      title: 'Settings',
      description: 'Setting the appearance of LayerView.',
      target: () => SettingsRef.current,
    },
    {
      title: 'Save',
      description: 'Save all your online-defined annotations as a json file.',
      target: () => SaveRef.current,
    },
    {
      title: 'Current Status',
      description: 'The Status of Current LayerView, including current Layer Species, Organs, Batch, Clustering, Embedding, Attributes, Batches, Spot number and Spot Radius',
      target: () => StatusRef.current,
    },
  ]

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

  const SpH5adLoader = (file) => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = (event) => {
        try {
          let h5info = H5adLoader(file, event, ['meta'])
          _setData(h5info.data)
          setTitle(meta?.st.dataset_id)
          setClusterOps(h5info.clusters)
          // console.log(h5info.clusters)
          if (h5info.clusters.map((item) => item.label).includes('batch'))
            setBatchName('batch')
          setEmbedOps(h5info.embdOps)
          resolve({
            data: h5info.data,
            title: h5info.title,
            clusters: h5info.clusters,
            embdOps: h5info.embdOps
          })
        } catch (error) {
          reject(error)
        }
      }
      reader.onloadend = () => {
        //console.log("Load sp data finished.")
      }
      reader.onerror = (error) => {
        reject(error)
      }
      reader.readAsArrayBuffer(file)
    })
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
    FileLoaderRef.current(file).then(() => {
      console.log("All data reloaded in LayerView.")
      toggleAnno("Upload")
    }).catch(error => {
      console.error('Error fetching blob in LayerView:', error)
    })
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

  useImperativeHandle(onRef, () => ({  // explode trigger for parent components
    "Trigger": toggleAnno, // Trigger for useEffect
    "Tour": setTourOpen, // Open the tutorial
    "Loading": setLoading, // set Loading status
    "Tip": setCurrTip, // set Loading Tips
  }))

  useEffect(() => {

    // console.log(commandRef.current)
    if (isInit) {
      var myChart = echarts.getInstanceByDom(chartRef.current)
      if (commandRef.current === "Reload") {
        setCompLoad((compLoad) => {
          let newCompLoad = { ...compLoad }
          newCompLoad['LayerView'] = true
          return newCompLoad
        })
        enterLoading(0, setLoading)
        setCurrTip(loadingTips[1])
        SpH5adLoader(spfile).then(() => {
          console.log("All data reloaded in LayerView.")
          toggleAnno("Upload")
        }).catch(error => {
          console.error('Error fetching blob in LayerView:', error)
        })
      }
      if (commandRef.current === "Upload") {
        enterLoading(0, setLoading)
        setCurrTip(loadingTips[2])
        // 1.set source
        let _dims = [...Object.keys(_data[0]), 'id']
        let source = _data.map((item, id) => {
          return [...Object.entries(item).map(([_, value]) => value), id]
        })
        symbolSizeRef.current = source.length > 5000 ? 4 : 6
        setItemSize(symbolSizeRef.current)
        setCellNum(source.length)
        // 2.set annotations and batches
        let defaultAnno = 'annotation'
        setClusterCur({ value: 0, label: defaultAnno, attr: 'categories' })
        prevCluster.current = { value: 0, label: defaultAnno, attr: 'categories' }
        let annotations = setItemGroup(source, _dims.indexOf(defaultAnno))
        // console.log(annotations)

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


        // 3.set embeddings
        let _embd = embedOps.map((item) => item.label).includes('spatial')
          ? embedOps.find((item) => item.label === 'spatial')
          : embedOps[0]
        let xName = `${_embd.label}_0`
        let yName = `${_embd.label}_1`
        let axis = setAxis(source, _dims, xName, yName)
        setEmbedCur(_embd)
        setEmbed3Dcur(_embd)

        // 4.set datasets
        let _datasets = setBatchDataset(source, _dims, defaultAnno, annotations, batches[0], _batchName)
        // console.log(_datasets)
        // 5.set 3D series
        let _series = []
        let _snum = 0
        let _3DSeries = {
          type: 'scatter3D',
          symbolSize: (symbolSizeRef.current / 2) > 2 ? (symbolSizeRef.current / 2) : 2,
          name: "batch",
          encode: {
            x: xName,
            y: yName,
            z: batchName,
            tooltip: [0],
          },
          itemStyle: {
            color: 'gray',
            opacity: 0.5,
          },
          datasetIndex: _snum,
          zlevel: -10,
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
            symbolSize: (symbolSizeRef.current / 2) > 2 ? (symbolSizeRef.current / 2) : 2,
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

        // console.log(axis)
        // console.log(_datasets)
        // console.log(_series)

        // 9.render charts
        myChart.setOption({
          title: [
            {
              text: title,
            },
          ],
          grid3D: {
            boxHeight: 135,
            top: '-8%',
            width: '60%',
            bottom: '20%',
          },
          grid: [{
            top: '14%',
            width: '40%',
            right: '1%',
            bottom: '18%',
          }],
          xAxis3D: axis.xAxis3D,
          yAxis3D: axis.yAxis3D,
          zAxis3D: { type: 'category', name: 'Batch', nameLocation: 'top' },
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
        quitLoading(0, setLoading)
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
        // console.log(_dataset)
        // set series
        let myVisualMap = []
        let _3DSeries = _series[0]
        _3DSeries.symbolSize = (itemSize / 2) > 2 ? itemSize : 2
        _3DSeries.itemStyle.opacity = itemOpacity / 2
        _3DSeries.encode = {
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
              symbolSize: (itemSize / 2) > 2 ? itemSize / 2 : 2,
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
            symbolSize: (itemSize / 2) > 2 ? (itemSize / 2) : 2,
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
            bottom: '0%',
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
                fontSize: token.fontSizeHeading4,
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
        //console.log(myChart.getOption())
      }
    } else {
      var myChart = echarts.init(chartRef.current) //init the echart container
      let axis = Axis.setEmptyAxis(0)
      let axis3D = Axis.setEmptyAxis3D(0)
      myChart.setOption({
        textStyle: {
          fontFamily: 'Arial'
        },
        tooltip: {},
        animation: false,
        title: [
          {
            text: title,
            top: '0%',
            left: 'center',
            textStyle: {
              fontSize: token.fontSizeHeading4,
            },
          },
        ],
        grid3D: {
          boxHeight: 150,
          top: '-5%',
          width: '60%',
          bottom: '25%',
        },
        xAxis3D: axis3D.xAxis3D,
        yAxis3D: axis3D.yAxis3D,
        zAxis3D: axis3D.zAxis3D,
        xAxis: axis.xAxis,
        yAxis: axis.yAxis,
        grid: [{
          top: '18%',
          width: '40%',
          right: '1%',
          bottom: '18%',
        }],
        legend: {
          type: 'scroll',
          animation: false,
          orient: 'horizontal',
          bottom: '0%',
          padding: 0,
          textStyle: {
            fontSize: 16,
          },
          pageIconColor: token.colorPrimaryActive,
          selector: [{ type: 'all', title: "All" }, { type: 'inverse', title: 'Inverse' }],
        },
        toolbox: {
          show: true,
          itemSize: 20,
          itemGap: 10,
          right: '0%',
          top: '0%',
          feature: {
            mark: { show: true },
            dataView: { show: true, readOnly: true },
            restore: { show: false },
            saveAsImage: {
              show: true,
              pixelRatio: 4,
              name: `LayerView_${title}`,

            },
            dataZoom: {},
          },
          iconStyle: {
            borderWidth: 1.5,
          },
        },
      })
      setInit(true)
    }
  }, [action])

  return (
    <Flex justify="center" gap='middle'>
      <Spin
        spinning={loading[0]}
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
          style={{ height: height, width: width, margin: margin, backgroundColor: token.colorBgContainer }}></div>
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
            <div ref={UploadRef}>
              <Dragger {...upLoadProps} style={{ marginTop: 20 }}>
                <p
                  className="ant-upload-drag-icon"
                  style={{ fontSize: 16, marginTop: 10, marginBottom: 10, alignItems: "center" }}>
                  <InboxOutlined />
                  <br />
                  Upload
                </p>
              </Dragger>
            </div>
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
                    <Button
                      icon={<ReloadOutlined />}
                      onClick={() => {
                        setItemSize(6)
                        setItemOpacity(0.8)
                        setyInv(false)
                        setxInv(false)
                        toggleAnno("Setting")
                      }}
                    >Reset</Button>
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
              <Button type="primary" block icon={<SettingOutlined />} ref={SettingsRef}>
                Settings
              </Button>
            </Popover>
            <Popover content={"Save the original file of current data"}>
              <Button type="primary" block
                icon={<CloudDownloadOutlined />}
                ref={SaveRef}
                onClick={() => {
                  let blob = spfile
                  let url = window.URL.createObjectURL(blob)
                  let a = document.createElement('a')
                  a.href = url
                  a.download = `sp_meta_${meta?.st.dataset_id}.h5ad`
                  document.body.appendChild(a)
                  a.click()

                  // 清理
                  window.URL.revokeObjectURL(url)
                  document.body.removeChild(a)
                }}>
                Save
              </Button>
            </Popover>
          </Space>
          <div ref={StatusRef}>
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
          </div>
          {/* <div>{JSON.stringify(meta)}</div> */}
        </Space>
      </ConfigProvider>
      <Tour open={tourOpen} onClose={() => setTourOpen(false)} steps={steps} />
    </Flex>
  )
}

SpScatter.defaultProps = {
  query: false,
  height: '35rem',
  width: '60rem',
  margin: '1rem',
}

SpScatter.propTypes = {
  spfile: PropTypes.object,
  setCompLoad: PropTypes.any,
  onRef: PropTypes.any,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
  meta: PropTypes.any,
  process: PropTypes.number,
}

export default SpScatter
