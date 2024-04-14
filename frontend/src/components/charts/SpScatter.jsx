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
import { CloudDownloadOutlined, FileAddOutlined, FileOutlined, InboxOutlined, SettingOutlined } from '@ant-design/icons'
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
  const [batchCur, setBatchCur] = useState({})
  const [embedCur, setEmbedCur] = useState('')
  const [clusterOps, setClusterOps] = useState({})
  const [embedOps, setEmbedOps] = useState({})
  const [batchOps, setBatchOps] = useState({})
  const [xInv, setxInv] = useState(false)
  const [yInv, setyInv] = useState(false)
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
  const setAxis = (source, dims, xName0, yName0,
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
        },
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

        // 2.set annotations and batches
        let defaultAnno = 'annotation'
        setClusterCur({ value: 0, label: defaultAnno, attr: 'categories' })
        prevCluster.current = { value: 0, label: defaultAnno, attr: 'categories' }
        let annotations = setItemGroup(source, _dims.indexOf(defaultAnno))

        let defaultBatch = 'batch'
        let batches = setItemGroup(source, _dims.indexOf(defaultBatch))
        setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
        setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))

        // 3.set embeddings
        let axis = setAxis(source, _dims, 'array_row', 'array_col')
        setEmbedCur('spatial')
        setEmbedOps([{ value: 0, label: 'spatial' }])
        let xName = 'array_row'
        let yName = 'array_col'

        // 4.set datasets
        let _datasets = setDataset(source, _dims, 'batch', batches)

        // 5.set 3D series
        let _series = []
        let _snum = 0
        let _3DSeries = {
          type: 'scatter3D',
          name: "batch",
          encode: {
            x: 'array_row',
            y: 'array_col',
            z: 'batch',
            tooltip: [0, 1, 2],
          },
          large: true,
          itemStyle: {
            color: 'gray',
            opacity: 0.3,
          },
          largeThreshold: 0,
          datasetIndex: _snum,
        }
        _series.push(_3DSeries)

        // 6.set batch series
        for (let bat of batches) {
          // seperate the _data into subgroups by cell-types
          _series.push({
            type: 'scatter',
            symbolSize: symbolSizeRef.current,
            xAxisIndex: 0,
            yAxisIndex: 0,
            name: bat,
            encode: {
              x: xName,
              y: yName,
              tooltip: [0, 1, 2, 3],
              itemName: bat,
            },
            emphasis: {
              focus: 'series',
            },
            large: true,
            largeThreshold: 5000,
            datasetIndex: _snum,
          })
          _snum = _snum + 1
        }
        console.log(_series)
        // 7.render charts
        myChart.setOption({
          tooltip: {},
          title: [
            {
              text: title,
              left: 'center',
              textStyle: {
                fontSize: 24,
              },
            },
          ],
          grid3D: {
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
          xAxis3D: {},
          yAxis3D: {},
          zAxis3D: { type: 'category' },
          xAxis: axis.xAxis,
          yAxis: axis.yAxis,
          grid: [{
            top: '15%',
            width: '40%',
            right: '1%',
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
          }],
          dataset: _datasets,
          series: _series,
        },
          {
            replaceMerge: ['dataset', 'series'],  // enable replaceMerge for datasets
          })
        console.log(myChart.getOption())
      }
      else if (commandRef.current === "Setting") {
        let option = myChart.getOption()
        let _source = option.dataset[0].source
        let _dims = option.dataset[0].dimensions
        let clusterIdx = _dims.indexOf(clusterCur.label)

        let _series = option.series
        // TODO: let batches change
      }
    } else {
      var myChart = echarts.init(chartRef.current, theme) //init the echart container

      let _dims = [...Object.keys(_data[0]), 'id']
      let source = _data.map((item, id) => {
        return [...Object.entries(item).map(([_, value]) => value), id]
      })
      symbolSizeRef.current = source.length > 5000 ? 2 : 4

      // set annotations
      let defaultAnno = 'leiden'
      setClusterCur({ value: 0, label: defaultAnno, attr: 'categories' })
      prevCluster.current = { value: 0, label: defaultAnno, attr: 'categories' }
      setClusterOps([{ value: 0, label: defaultAnno, attr: 'categories' }])
      let annotations = setItemGroup(source, _dims.indexOf('leiden'))


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
          large: true,
          itemStyle: {
            color: 'gray',
            opacity: 0.3,
          },
          largeThreshold: 0,
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
            left: 'center',
            textStyle: {
              fontSize: 24,
            },
          },
        ],
        grid3D: {
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
        xAxis3D: {},
        yAxis3D: {},
        zAxis3D: {},
        xAxis: axis.xAxis,
        yAxis: axis.yAxis,
        grid: [{
          top: '15%',
          width: '40%',
          right: '1%',
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
        }],
        dataset: _datasets,
        series: _series,
      })
      setInit(true)
    }
  }, [theme, action])

  return (
    <Flex justify="center" gap='middle'>
      <div
        ref={chartRef}
        className="chart"
        //the target DOM container needs height and width
        style={{ height: height, width: width, margin: margin }}></div>
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
        <Button type="primary" block icon={<SettingOutlined />}>
          Settings
        </Button>
        <Button type="primary" block icon={<CloudDownloadOutlined />}>
          Save
        </Button>
        <Form
          name="Settings"
          size="middle"
          labelCol={{
            span: 6,
          }}
          wrapperCol={{
            span: 17,
          }}>
          {/* <Form.Item
            label="Batch"
            labelCol={{
              span: 7,
            }}>
            {' '}
            <Select
              labelInValue
              placeholder="Batch"
              style={{
                width: '100%',
              }}
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
          </Form.Item> */}
        </Form>
      </Space>
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
