import React, { useRef, useState, useEffect } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import vega_20 from '../theme/vega_20'
import {
  GraphicComponent,
  GridComponent,
  LegendComponent,
  TooltipComponent,
} from 'echarts/components'
import { LineChart } from 'echarts/charts'
import { UniversalTransition } from 'echarts/features'
import { CanvasRenderer } from 'echarts/renderers'
import axios from 'axios'
import saveAs from 'file-saver'
import H5adLoader from '../utils/H5adLoader'
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
  theme,
  Col,
  Row,
  Form,
  Flex,
  Spin
} from 'antd'
import {
  CheckOutlined,
  CloseOutlined,
  CloudDownloadOutlined,
  DeleteOutlined,
  FormOutlined,
  ReloadOutlined,
  SettingOutlined,
  SlidersOutlined,
  UploadOutlined,
  InboxOutlined,
} from '@ant-design/icons'
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

const PairView = ({ height, width, margin }) => {
  const [api, contextHolder] = notification.useNotification()
  const [isInit, setInit] = useState(false) // whether echart object is inited
  const chartRef = useRef(null) // current DOM container
  const [_scdata, _setScData] = useState(
    require('../../assets/data/CID4971-umap-sc.json')
  )
  const [_spdata, _setSpData] = useState(
    {}
  )
  const [seriesArray, setSeriesArray] = useState([])
  const [nameArray, setNameArray] = useState([])
  const snumRef = useRef()
  const datasetRef = useRef()
  const seriesRef = useRef()
  const commandRef = useRef()
  const nameRef = useRef()
  const brushRef = useRef([])
  const [brushArray, setBrushArray] = useState([])
  const [inputValue, setInputValue] = useState('')
  const [allowConfirm, setAllowConfirm] = useState(false)
  const [action, setAction] = useState(0)
  const [deleteValue, setDeleteValue] = useState([])
  const [deleteOpen, setDeleteOpen] = useState(false)
  const [refineValue, setRefineValue] = useState([])
  const [refineOpen, setRefineOpen] = useState(false)
  const itemGroupRef = useRef([])
  const [title, setTitle] = useState('CID4971-umap-sc')
  const FileLoaderRef = useRef('')
  const symbolSizeRef = useRef('')
  const [loadings, setLoadings] = useState([])
  const [Uploading, setUploading] = useState(false)
  const brushMode = useRef('Select')
  const [brushModeState, setBrushModeState] = useState('Select')

  // Item size and Opacity for Single-cell and Spatial
  const [itemSizeSc, setItemSizeSc] = useState(2)
  const [itemOpacitySc, setItemOpacitySc] = useState(0.8)
  const [itemSizeSp, setItemSizeSp] = useState(4)
  const [itemOpacitySp, setItemOpacitySp] = useState(0.8)
  const prevCluster = useRef('')

  // Clusters and Embedings for Single-cell
  const [clusterCurSc, setClusterCurSc] = useState({})
  const [embedCurSc, setEmbedCurSc] = useState('')
  const [clusterOpsSc, setClusterOpsSc] = useState({})
  const [embedOpsSc, setEmbedOpsSc] = useState({})
  const [annoCurSc, setAnnoCurSc] = useState([])

  // Clusters and Embedings for ST
  const [clusterCurSp, setClusterCurSp] = useState({})
  const [embedCurSp, setEmbedCurSp] = useState('')
  const [clusterOpsSp, setClusterOpsSp] = useState({})
  const [embedOpsSp, setEmbedOpsSp] = useState({})
  const [annoCurSp, setAnnoCurSp] = useState([])
  const stAnnoLength = useRef(0)

  const [xInv, setxInv] = useState(false)
  const [yInv, setyInv] = useState(false)
  const [cellNum, setCellNum] = useState(0)
  const [refineOption, setRefineOption] = useState([
    { value: 0, label: 'Eager Refiner' },
    { value: 1, label: 'Lazy Refiner' },
    { value: 2, label: 'LabelPropagation' },
    { value: 3, label: 'LabelSpreading' },
  ])
  const { token } = useToken()

  const setItemGroup = (source, annoIdx, type = 'categories', doSet = true) => {
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
    if (doSet) {
      itemGroupRef.current = itemGroup
    }
    return annotations
  }

  const setBrushedMap = (option, brushed, doPop = true, mode) => {
    let __datasets = option.dataset
    let __series = option.series
    let __encode = __series[0].encode
    let brushedNum = option.series.length + 1
    let popDataset = {}
    if (doPop) {
      popDataset = __datasets.pop()
    } else {
      // default name for first brush
      nameRef.current = `annotation ${brushedNum - stAnnoLength.current - itemGroupRef.current.length}`
      seriesRef.current = {
        type: 'scatter',
        xAxisIndex: 0,
        yAxisIndex: 0,
        symbolSize: 10,
        name: nameRef.current,
        encode: {
          x: __encode.x,
          y: __encode.y,
        },
        large: true,
        largeThreshold: 5000,
        datasetIndex: option.dataset.length,
      }
      __series.push(seriesRef.current)
    }
    let newBrushed = brushed
    if (mode === 'Erase') {
      let prevBrushed = popDataset.transform.config.or.map((item) => {
        return item['=']
      })
      newBrushed = prevBrushed.filter((item) => !brushed.includes(item))
    } else if (mode === 'Draw' && doPop) {
      let prevBrushed = popDataset.transform.config.or.map((item) => {
        return item['=']
      })
      newBrushed = [...new Set([...prevBrushed, ...brushed])]
    } else {
      newBrushed = brushed
    }
    brushRef.current = newBrushed
    datasetRef.current = {
      transform: {
        type: 'filter',
        config: {
          or: newBrushed.map((item) => ({
            dimension: 'id',
            '=': item,
          })),
        },
      },
      fromDatasetIndex: 0,
    }
    __datasets.push(datasetRef.current)
    if (doPop) {
      return {
        dataset: __datasets,
      }
    } else {
      return {
        dataset: __datasets,
        series: __series,
      }
    }
  }

  const setEmptyAxis = (gridIndex) => ({
    xAxis: {
      id: gridIndex,
      gridIndex: gridIndex,
      name: gridIndex,
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
    },
    yAxis:
    {
      id: gridIndex,
      gridIndex: gridIndex,
      name: gridIndex,
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
    },
  })

  const setAxis = (
    source,
    dims,
    xName,
    yName,
    gridIndex,
  ) => {
    // the axis of left chart
    let xIdx = dims.indexOf(xName)
    let yIdx = dims.indexOf(yName)
    let embd_x = source.map((item) => {
      return item[xIdx]
    })
    let embd_y = source.map((item) => {
      return item[yIdx]
    })
    return {
      xAxis: {
        id: gridIndex,
        gridIndex: gridIndex,
        name: xName,
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
        min: Math.ceil(Math.min.apply(null, embd_x) - 1.5),
        max: Math.ceil(Math.max.apply(null, embd_x) + 1.5),
      },
      yAxis:
      {
        id: gridIndex,
        gridIndex: gridIndex,
        name: yName,
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
        min: Math.ceil(Math.min.apply(null, embd_y) - 1.5),
        max: Math.ceil(Math.max.apply(null, embd_y) + 1.5),
      },
    }
  }

  const SpH5adLoader = (file) => {
    const reader = new FileReader()
    reader.onload = (event) => {
      let h5info = H5adLoader(file, event)
      _setSpData(h5info.data)
      setTitle(h5info.title)
      setClusterOpsSp(h5info.clusters)
      setEmbedOpsSp(h5info.embdOps)
    }
    reader.onloadend = () => {
      toggleAnno('Upload')
    }
    reader.readAsArrayBuffer(file)
  }

  const setDataset = (source, dims, annoName, annotations, fromIndex = 0) => {
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
        fromDatasetIndex: fromIndex,
      })
    } else {
      for (let anno of annotations) {
        _datasets.push({
          // 这个 dataset 的 index 是 `1`。
          transform: {
            type: 'filter',
            config: { dimension: annoName, value: anno },
          },
          fromDatasetIndex: fromIndex,
        })
      }
    }
    return _datasets
  }

  const enterLoading = (index) => {
    setLoadings((prevLoadings) => {
      let newLoadings = [...prevLoadings]
      newLoadings[index] = true
      return newLoadings
    })
    setTimeout(() => {
      setLoadings((prevLoadings) => {
        let newLoadings = [...prevLoadings]
        newLoadings[index] = false
        return newLoadings
      })
    }, 60000)
  }

  const quitLoading = (index) => {
    setLoadings((prevLoadings) => {
      let newLoadings = [...prevLoadings]
      newLoadings[index] = false
      return newLoadings
    })
  }

  useEffect(() => {
    if (isInit) {
      // get the echart container
      var myChart = echarts.getInstanceByDom(chartRef.current)
      let _series = seriesArray
      if (commandRef.current === 'Upload') {
        let _scdims = [...Object.keys(_scdata[0]), 'id']
        let _scsource = _scdata.map((item, id) => {
          return [...Object.entries(item).map(([_, value]) => value), id]
        })
        let _spdims = [...Object.keys(_spdata[0]), 'id']
        let _spsource = _spdata.map((item, id) => {
          return [...Object.entries(item).map(([_, value]) => value), id]
        })
        symbolSizeRef.current = _scsource.length > 5000 ? 2 : 4
        setCellNum(_spsource.length)

        // set sc annotations
        let scdefaultAnno
        if (_scdims.includes('annotation')) {
          scdefaultAnno = 'annotation'
        } else {
          scdefaultAnno = _scdims[0]
        }
        setClusterCurSc({ value: 0, label: scdefaultAnno, attr: 'categories' })
        prevCluster.current = { value: 0, label: scdefaultAnno, attr: 'categories' }
        setClusterOpsSc([{ value: 0, label: scdefaultAnno, attr: 'categories' }, { value: 1, label: 'umap_0', attr: 'numeric' }])
        let _scannotations = setItemGroup(_scsource, _scdims.indexOf(scdefaultAnno), 'categories', true)
        setAnnoCurSc(_scannotations)

        // set sp annotations
        let spdefaultAnno
        if (_spdims.includes('annotation')) {
          spdefaultAnno = 'annotation'
        } else {
          spdefaultAnno = 'leiden'
        }
        setClusterCurSp({ value: 0, label: spdefaultAnno, attr: 'categories' })
        let _spannotations = setItemGroup(_spsource, _spdims.indexOf(spdefaultAnno), 'categories', false)
        setAnnoCurSp(_spannotations)
        stAnnoLength.current = _spannotations.length

        // set sc embeddings
        let _scaxis = setAxis(_scsource, _scdims, 'umap_0', 'umap_1', 0)
        setEmbedCurSc('umap')
        setEmbedOpsSc([{ value: 0, label: 'umap' }])

        // set sp embeddings
        let _spaxis = setAxis(_spsource, _spdims, 'array_row', 'array_col', 1)
        setEmbedCurSp('spatial')


        // set datasets
        let _scdatasets = setDataset(_scsource, _scdims, scdefaultAnno, _scannotations, 0)
        let _spdatasets = setDataset(_spsource, _spdims, spdefaultAnno, _spannotations, _scdatasets.length)
        let _datasets = [..._scdatasets, ..._spdatasets]
        let _series = []
        let _snum = 0

        // set sc annotations
        for (let anno of _scannotations) {
          _snum = _snum + 1
          // seperate the _data into subgroups by cell-types
          _series.push({
            type: 'scatter',
            symbolSize: symbolSizeRef.current,
            xAxisIndex: 0,
            yAxisIndex: 0,
            name: anno,
            encode: {
              x: 'umap_0',
              y: 'umap_1',
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
        // set sp annotations
        _snum = _snum + 1
        for (let anno of _spannotations) {
          _snum = _snum + 1
          // seperate the _data into subgroups by cell-types
          _series.push({
            type: 'scatter',
            xAxisIndex: 1,
            yAxisIndex: 1,
            symbolSize: symbolSizeRef.current + 2,
            name: anno,
            encode: {
              x: 'array_row',
              y: 'array_col',
              itemName: anno,
            },
            emphasis: {
              focus: 'series',
            },
            datasetIndex: _snum,
          })
        }
        snumRef.current = _snum

        // set visualMap
        myChart.setOption({
          color: vega_20,
          title: [
            {
              text: title,
              left: 'center',
              textStyle: {
                fontSize: 24,
              },
            },
          ],
          xAxis: [_scaxis.xAxis, _spaxis.xAxis],
          yAxis: [_scaxis.yAxis, _spaxis.yAxis],
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
          dataset: _datasets,
          series: _series,
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
          },
        }, {
          replaceMerge: ['dataset', 'series', 'visualMap'],  // enable replaceMerge for datasets
        })
      }
      if (commandRef.current === 'scConfigs') {
        let option = myChart.getOption()

        // set sc source
        let _scsource = option.dataset[0].source
        let _scdims = option.dataset[0].dimensions
        let scclusterIdx = _scdims.indexOf(clusterCurSc.label)

        // set scaxis
        let xName0 = embedCurSc + '_0'
        let yName0 = embedCurSc + '_1'
        let scaxis = setAxis(_scsource, _scdims, xName0, yName0, 0)

        // set datasets if cluster changes
        // datasets
        // |-|-------|-|------------|----|
        // |--sclen--| |---stAnno---|
        // |-----prevAnnolen--------|

        // series
        // |--------|------------|
        // |-sclen-1|---stAnno---|
        let _sclen = annoCurSc.length + 1
        if (annoCurSc.length === 0) {
          _sclen = _sclen + 1
        }
        let clu = clusterCurSc.label
        let annotations = setItemGroup(_scsource, scclusterIdx, clusterCurSc.attr, true)
        let prevAnnoLen = stAnnoLength.current + _sclen + 1
        if(stAnnoLength.current === 0){
          prevAnnoLen = prevAnnoLen + 1
        }
        let annoLenOffset = _sclen - 1 - annotations.length
        setAnnoCurSc(annotations)
        if (annotations.length === 0) {
          annoLenOffset = annoLenOffset - 1
        }
        let _scdatasets = setDataset(_scsource, _scdims, clu, annotations, 0)
        let _spdatasets = option.dataset.slice(_sclen, prevAnnoLen)
        for (let i = 1; i < _spdatasets.length; i++) {
          _spdatasets[i].fromDatasetIndex = _scdatasets.length
        }
        let _datasets = [..._scdatasets, ..._spdatasets, ...option.dataset.slice(prevAnnoLen)]

        // set scseries
        let _snum = 0
        let _scseries = []
        let myVisualMap = []
        if (clusterCurSc.attr === 'categories') {
          for (let anno of annotations) {
            // seperate the _data into subgroups by cell-types
            _snum = _snum + 1
            _scseries.push({
              type: 'scatter',
              xAxisIndex: 0,
              yAxisIndex: 0,
              symbolSize: itemSizeSc,
              name: anno,
              encode: {
                // annotations are displayed in left chart
                x: xName0,
                y: yName0,
                itemName: anno,
              },
              emphasis: {
                focus: 'series',
              },
              itemStyle: {
                opacity: itemOpacitySc,
              },
              datasetIndex: _snum,
            })
          }
        } else {
          _snum = _snum + 1
          _scseries.push({
            type: 'scatter',
            name: clusterCurSc.label,
            symbolSize: itemSizeSc,
            encode: {
              // annotations are displayed in left chart
              x: xName0,
              y: yName0,
              tooltip: [
                0,
                _scdims.indexOf(xName0),
                _scdims.indexOf(yName0),
                scclusterIdx,
                _scdims[_scdims.length - 1],
              ],
            },
            itemStyle: {
              opacity: itemOpacitySc,
            },
            datasetIndex: _snum,
          })
          myVisualMap.push({
            show: true,
            id: clu + '_sc',
            calculable: true,
            dimension: clu,
            seriesIndex: 0,
            left: '10%',
            top: '5%',
            orient: 'horizontal',
            precision: 2,
            calculable: true,
            min: Math.min.apply(
              null,
              _scsource.map((item) => item[scclusterIdx])
            ),
            max: Math.max.apply(
              null,
              _scsource.map((item) => item[scclusterIdx])
            ),
            inRange: {
              color: ['#808080', '#FF3300'],
              opacity: [0.5, 1],
            },
            text: ['', clusterCurSc.label],
            textGap: 20,
            textStyle: {
              fontSize: 16,
            },
          })
        }
        // append brush clusters and sp clusters
        _snum = _snum + 1
        let manualSeries = option.series.slice(_sclen - 1)
        for (let i = 0; i < manualSeries.length; i++) {
          _snum = _snum + 1
          manualSeries[i].datasetIndex = _snum
        }
        let _series = [..._scseries, ...manualSeries]
        snumRef.current = _datasets.length - 1 // the sum of default and manual annotation

        // set VisualMap
        let preVisualMap = option.visualMap
        let _visualMap
        if (option.visualMap.length === 0) {
          _visualMap = myVisualMap
        } else if (option.visualMap.length === 2) {
          preVisualMap[1].seriesIndex = _scseries.length
          _visualMap = [...myVisualMap, preVisualMap[1]]
        } else if (option.visualMap.length === 1) { //check id is _sc or _sp
          let _strs = preVisualMap[0].id.split("_")
          let id = _strs[_strs.length - 1]
          if (id === 'sp') {
            preVisualMap[0].seriesIndex = _scseries.length
            _visualMap = [...myVisualMap, ...preVisualMap]
          } else {
            _visualMap = myVisualMap
          }
        }
        console.log("sc", _visualMap)
        myChart.setOption(
          {
            xAxis: [scaxis.xAxis, option.xAxis[1]],
            yAxis: [scaxis.yAxis, option.yAxis[1]],
            dataset: _datasets,
            series: _series,
            visualMap: _visualMap,
          },
          {
            replaceMerge: ['dataset', 'series', 'visualMap'],
          }
        )
        prevCluster.current = clusterCurSc
        console.log(myChart.getOption())
      }
      if (commandRef.current === 'spConfigs') {
        let option = myChart.getOption()

        // set sp source
        let _sclen = annoCurSc.length + 1
        if (annoCurSc.length === 0) {
          _sclen = _sclen + 1
        }
        let _spsource = option.dataset[_sclen].source
        let _spdims = option.dataset[_sclen].dimensions

        // set sp Axis
        let xName = embedCurSp + '_0'
        let yName = embedCurSp + '_1'
        if(embedCurSp === 'spatial' && _spdims.includes("array_row") && _spdims.includes("array_col")){
          xName = 'array_row'
          yName = 'array_col'
        }
        let spaxis = setAxis(_spsource, _spdims, xName, yName, 1)

        // set datasets
        // |-|-------|-|------------|----|
        // |--sclen--| |---stAnno---|
        // |-----prevAnnolen--------|

        // series
        // |--------|------------|
        // |-sclen-1|---stAnno---|
        let clu = clusterCurSp.label
        let cluIdx = _spdims.indexOf(clusterCurSp.label)
        let prevAnnoLen = stAnnoLength.current + _sclen + 1
        if(stAnnoLength.current === 0){
          prevAnnoLen = prevAnnoLen + 1
        }
        let annotations = setItemGroup(_spsource, cluIdx, clusterCurSp.attr, false)  // spconfig do not change the item group
        let annoLenOffset = stAnnoLength.current - annotations.length
        if(stAnnoLength.current === 0){
          annoLenOffset = annoLenOffset + 1
        }
        stAnnoLength.current = annotations.length
        setAnnoCurSp(annotations)
        let _spdatasets = setDataset(_spsource, _spdims, clu, annotations, _sclen)
        let _datasets = [...option.dataset.slice(0, _sclen), ..._spdatasets, ...option.dataset.slice(prevAnnoLen)]
        console.log("dataset-prevAnnolen",option.dataset.slice(prevAnnoLen))
        // set series and VisualMap
        let _spSeries = []
        let _snum = _sclen
        let myVisualMap = []
        if (clusterCurSp.attr === 'categories') {
          for (let anno in annotations) {
            _snum = _snum + 1
            _spSeries.push({
              type: 'scatter',
              symbolSize: itemSizeSp,
              xAxisIndex: 1,
              yAxisIndex: 1,
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
                opacity: itemOpacitySp
              },
              datasetIndex: _snum,
            })
          }
        } else {
          annoLenOffset = annoLenOffset - 1  // for numeric, annotations.length = 0, but it in fact uses a series num, so the offset should minus 1.
          _spSeries.push({
            type: 'scatter',
            name: clu,
            xAxisIndex: 1,
            yAxisIndex: 1,
            symbolSize: itemSizeSp,
            encode: {
              // annotations are displayed in left chart
              x: xName,
              y: yName,
            },
            itemStyle: {
              opacity: itemOpacitySp,
            },
            datasetIndex: _snum,
          })

          myVisualMap.push({
            show: true,
            id: clu + '_sp',
            calculable: true,
            dimension: clu,
            seriesIndex: _sclen - 1,
            right: '0%',
            top: '5%',
            orient: 'horizontal',
            precision: 2,
            calculable: true,
            min: Math.min.apply(
              null,
              _spsource.map((item) => item[cluIdx])
            ),
            max: Math.max.apply(
              null,
              _spsource.map((item) => item[cluIdx])
            ),
            inRange: {
              color: ['#808080', '#FF3300'],
              opacity: [0.8, 1],
            },
            text: [clu, ''],
            textGap: 20,
            textStyle: {
              fontSize: 16,
            },
          })
        }
        // append annoated clusters
        let manualSeries = option.series.slice(prevAnnoLen - 2)
        for (let i = 0; i < manualSeries.length; i++) {
          manualSeries[i].datasetIndex =
            manualSeries[i].datasetIndex - annoLenOffset
        }
        console.log(manualSeries)
        let _series = [...option.series.slice(0, _sclen-1), ..._spSeries, ...manualSeries]

        // set VisualMap
        let preVisualMap = option.visualMap
        let _visualMap
        if (option.visualMap.length === 0) {
          _visualMap = myVisualMap
        } else if (option.visualMap.length === 2) {
          _visualMap = [preVisualMap[0], ...myVisualMap]
        } else if (option.visualMap.length === 1) { //check id is _sc or _sp
          let _strs = preVisualMap[0].id.split("_")
          let id = _strs[_strs.length - 1]
          if (id === 'sc') {
            _visualMap = [...preVisualMap, ...myVisualMap]
          } else {
            _visualMap = myVisualMap
          }
        }
        console.log("sp", _visualMap)

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
          xAxis: [option.xAxis[0], spaxis.xAxis],
          yAxis: [option.yAxis[0], spaxis.yAxis],
          series: _series,
          dataset: _datasets,
          visualMap: _visualMap,
        },
          {
            replaceMerge: ['series', 'dataset', 'visualMap'],
          })
        console.log(myChart.getOption())
      }
      if (commandRef.current === 'Confirm') {
        snumRef.current = snumRef.current + 1
        console.log(`seriesNum ${snumRef.current}`)
        _series.push(seriesRef.current)
        setSeriesArray(_series)
        setBrushArray([
          ...brushArray,
          {
            name: nameRef.current,
            data: brushRef.current,
            value: snumRef.current - itemGroupRef.current.length,
          },
        ])
        setNameArray([
          ...nameArray,
          {
            value: snumRef.current - itemGroupRef.current.length,
            label: nameRef.current,
          },
        ])
        setAllowConfirm(false)
        api['success']({
          message: 'Annotation Confirmed',
          description: `Your Annotation \'${nameRef.current}\' is confirmed.`,
          placement: 'topRight',
        })
        datasetRef.current = null
        seriesRef.current = null
        commandRef.current = null
        brushRef.current = []
        setBrushModeState('Select')
        brushMode.current = 'Select'
      }
      if (commandRef.current === 'Rename') {
        seriesRef.current.name = inputValue
        let prevName = nameRef.current
        nameRef.current = inputValue
        _series.push(seriesRef.current)
        myChart.setOption({
          series: _series,
        })
        setInputValue('')
        api.info({
          message: 'Annotation Renamed',
          description: `Previous Annotation \'${prevName}\' is changed to \'${nameRef.current}\'.`,
          placement: 'topRight',
        })
        commandRef.current = null
      }
      if (commandRef.current === 'Refine') {
        let starttime = Date.now()
        axios
          .post('http://localhost:5522/refine', {
            data: {
              name: `resources/${title}.h5ad`,
              anno: brushRef.current,
              refiner: refineValue.value,
              starttime: starttime,
            },
            headers: {
              'Content-Type': 'application/json',
            },
          })
          .then((response) => {
            if (response.data.success) {
              let endtime = response.data.endtime
              let refScatter = response.data.refined
              let option = myChart.getOption()
              let newOption = setBrushedMap(option, refScatter, true, 'Select')
              myChart.setOption(newOption)
              api.success({
                message: `Annotation Refined in ${(endtime - starttime) / 1000}s`,
                description: response.data.message,
                placement: 'topRight',
              })
            } else {
              let endtime = response.data.endtime
              api.warning({
                message: `Annotation does not Refined in ${(endtime - starttime) / 1000}s`,
                description: response.data.message,
                placement: 'topRight',
              })
            }
            setLoadings((prevLoadings) => {
              let newLoadings = [...prevLoadings]
              newLoadings[0] = false
              return newLoadings
            })
          }).catch((error) => {
            api.error({
              message: `Annotation failed for ${error}`,
              placement: 'topRight',
            })
          })
        quitLoading(0)
      }
      if (commandRef.current === 'Delete') {
        let option = myChart.getOption()
        let __series = option.series
        setNameArray(
          nameArray.filter((item) => !deleteValue.includes(item.value))
        )
        setBrushArray(
          brushArray.filter((item) => !deleteValue.includes(item.value))
        )
        let _indices = deleteValue.map(
          (item) => item + itemGroupRef.current.length
        )
        __series = __series.filter(
          (item) => !_indices.includes(item.datasetIndex)
        )
        myChart.setOption(
          {
            series: __series,
          },
          {
            replaceMerge: ['series'],
          }
        )
        setSeriesArray(__series)
        setDeleteValue([])
        console.log(myChart.getOption())
        api.warning({
          message: 'Annotation Deleted',
          description:
            'Previous annotation is deleted, but you can recover them until leaving this page.',
          placement: 'topRight',
        })
      }
      if (commandRef.current === 'Download') {
        let jsonFile = {}
        jsonFile.columns = ['barcode', 'umap_0', 'umap_1', 'annotation']
        jsonFile.data = []
        for (let item of brushArray) {
          let itemData = item.data.map((index) => [
            _scdata[index]['index'],
            _scdata[index]['umap_0'],
            _scdata[index]['umap_1'],
            item.name,
          ])
          jsonFile.data = [...jsonFile.data, ...itemData]
        }
        const fileStream = JSON.stringify(jsonFile)
        const blob = new Blob([fileStream])
        saveAs(blob, `Annotation-${new Date().getTime()}.json`)
      }
    } else {
      // init the echart container
      var myChart = echarts.init(chartRef.current)
      axios
        .get('http://localhost:5522/query', {
          responseType: 'blob',
        })
        .then((response) => {
          let blob = response.data
          const file = new File([blob], "sp1_meta.h5ad")
          SpH5adLoader(file)
        })
        .catch(error => {
          console.error('Error fetching blob:', error);
        })
      let _scaxis = setEmptyAxis(0)
      let _spaxis = setEmptyAxis(1)
      myChart.setOption({
        tooltip: {},
        brush: {
          toolbox: ['rect', 'polygon', 'clear'],
          outOfBrush: {
            colorAlpha: 0.3,
          },
          xAxisIndex: 0,
          throttleType: 'debounce',
          throttleDelay: 1000,
        },
        toolbox: {
          show: true,
          itemSize: 20,
          itemGap: 10,
          feature: {
            bursh: {
              type: ['rect', 'polygon', 'clear'],
            },
            mark: { show: true },
            dataView: { show: true, readOnly: true },
            restore: { show: true },
            saveAsImage: { show: true },
            dataZoom: {},
          },
          iconStyle: {
            borderWidth: 1.5,
          },
          top: '0%',
          left: '0%',
        },
        xAxis: [_scaxis.xAxis, _spaxis.xAxis],
        yAxis: [_scaxis.yAxis, _spaxis.yAxis],
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
      })
      setInit(true)

      myChart.on('brushSelected', function (params) {
        let brushed = []
        let brushComponent = params.batch[0]
        let selectedIdx = brushComponent.selected.map((item) => {
          return item.seriesIndex
        })
        let selectedInd = brushComponent.selected.map((item) => {
          return item.dataIndex
        })
        for (let i = 0; i < selectedIdx.length; i++) {
          // search for each seriesIndex
          let sIdx = selectedIdx[i]
          let dIdx = selectedInd[i]
          for (let j = 0; j < dIdx.length; j++) {
            let dataIndex = dIdx[j]
            brushed.push(itemGroupRef.current[sIdx][dataIndex])
          }
        }

        if (brushed.length > 0) {
          // trigger the brush command in useEffect
          let option = myChart.getOption()
          let newOption = {}
          if (brushRef.current.length > 0) {
            newOption = setBrushedMap(option, brushed, true, brushMode.current)
          } else {
            newOption = setBrushedMap(option, brushed, false, brushMode.current)
          }
          myChart.setOption(newOption)
          setAllowConfirm(true)
        }
      })
    }
  }, [action])

  const toggleAnno = (command) => {
    commandRef.current = command
    setAction(action + 1)
  }

  const ScH5adLoader = (file) => {
    const reader = new FileReader()
    reader.onload = (event) => {
      let h5info = H5adLoader(file, event)
      _setScData(h5info.data)
      setTitle(h5info.title)
      setClusterOpsSc(h5info.clusters)
      setEmbedOpsSc(h5info.embdOps)
    }
    reader.onloadend = () => {
      toggleAnno('Upload')
      setUploading(false)
    }
    reader.readAsArrayBuffer(file)
  }

  const JsonLoader = (file) => {
    let newTitle = file.name.replace(/\.json$/, '')
    setTitle(newTitle)
    const reader = new FileReader()
    reader.onload = (event) => {
      const jsonData = JSON.parse(event.target.result)
      _setScData(jsonData)
      toggleAnno('Upload')
    }
    reader.readAsText(file)
  }

  const CsvLoader = (file) => {
    let newTitle = file.name.replace(/\.json$/, '')
    setTitle(newTitle)
    const reader = new FileReader()
    reader.onload = (event) => {
      const jsonData = JSON.parse(event.target.result)
      _setScData(jsonData)
      toggleAnno('Upload')
    }
    reader.readAsText(file)
  }

  return (
    <div>
      {contextHolder}

      <Flex justify="center" gap='middle'>
        <Spin spinning={Uploading} size="large">
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
          <Space direction="vertical" size="small">
            <Space direction="vertical" size="small">
              <Popover
                title={<>Setting Figure Options for single-cell view(left):</>}
                placement="topLeft"
                style={{ width: 800 }}
                content={
                  <Space size="small" direction="vertical">
                    <Form
                      name="Settings"
                      size="large"
                      labelCol={{
                        span: 6,
                      }}
                      wrapperCol={{
                        span: 17,
                      }}>
                      <Form.Item label="ItemSize">
                        {' '}
                        <Slider
                          min={0}
                          max={10}
                          marks={{ 0: '0', 5: '5', 10: '10' }}
                          onChange={(value) => {
                            setItemSizeSc(value)
                          }}
                          value={typeof itemSizeSc === 'number' ? itemSizeSc : 2}
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
                            setItemOpacitySc(value)
                          }}
                          value={
                            typeof itemOpacitySc === 'number' ? itemOpacitySc : 0.8
                          }
                        />
                      </Form.Item>
                      <Form.Item
                        label="Label Loc"
                        labelCol={{
                          span: 7,
                        }}>
                        {' '}
                        <Segmented
                          options={[
                            { label: 'on Data', value: 'on Data' },
                            { label: 'Bottom', value: 'Bottom' },
                            { label: 'Hide', value: 'Hide' },
                          ]}
                        />
                      </Form.Item>
                      <Form.Item
                        label="Brush Mode"
                        labelCol={{
                          span: 8,
                        }}>
                        {' '}
                        <Segmented
                          options={[
                            { label: 'Select', value: 'Select' },
                            { label: 'Draw', value: 'Draw' },
                            { label: 'Erase', value: 'Erase' },
                          ]}
                          defaultValue="Select"
                          onChange={(value) => {
                            brushMode.current = value
                            setBrushModeState(value)
                          }}
                          value={brushModeState}
                        />
                      </Form.Item>
                      <Form.Item
                        label="Cluster"
                        labelCol={{
                          span: 7,
                        }}>
                        {' '}
                        <Select
                          labelInValue
                          placeholder="Cluster"
                          style={{
                            width: '100%',
                          }}
                          placement="topLeft"
                          options={clusterOpsSc}
                          value={clusterCurSc}
                          onChange={(target) => {
                            setClusterCurSc({
                              value: target.value,
                              label: target.label,
                              attr: clusterOpsSc[target.value].attr,
                            })
                            toggleAnno('scConfigs')
                          }}
                        />
                      </Form.Item>
                      <Form.Item
                        label="Embedding"
                        labelCol={{
                          span: 9,
                        }}>
                        {' '}
                        <Select
                          labelInValue
                          placeholder="Embedding"
                          style={{
                            width: '100%',
                          }}
                          placement="topLeft"
                          options={embedOpsSc}
                          value={embedCurSc}
                          onChange={(value) => {
                            setEmbedCurSc(value.label)
                            toggleAnno('scConfigs')
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
                          toggleAnno('scConfigs')
                        }}>
                        Apply
                      </Button>
                    </Space>
                  </Space>
                }
                trigger="click">
                <Button type="primary" block icon={<SettingOutlined />}>
                  scConfigs
                </Button>
              </Popover>
              <Popover
                title={<>Options for Spatial View(right):</>}
                placement="topLeft"
                style={{ width: 800 }}
                content={
                  <Space size="small" direction="vertical">
                    <Form
                      name="Settings"
                      size="large"
                      labelCol={{
                        span: 6,
                      }}
                      wrapperCol={{
                        span: 17,
                      }}>
                      <Form.Item label="ItemSize">
                        {' '}
                        <Slider
                          min={0}
                          max={10}
                          marks={{ 0: '0', 5: '5', 10: '10' }}
                          onChange={(value) => {
                            setItemSizeSp(value)
                          }}
                          value={typeof itemSizeSp === 'number' ? itemSizeSp : 2}
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
                            setItemOpacitySp(value)
                          }}
                          value={
                            typeof itemOpacitySp === 'number' ? itemOpacitySp : 0.8
                          }
                        />
                      </Form.Item>
                      <Form.Item
                        label="Label Loc"
                        labelCol={{
                          span: 7,
                        }}>
                        {' '}
                        <Segmented
                          options={[
                            { label: 'on Data', value: 'on Data' },
                            { label: 'Bottom', value: 'Bottom' },
                            { label: 'Hide', value: 'Hide' },
                          ]}
                        />
                      </Form.Item>
                      <Form.Item
                        label="Inverse"
                        labelCol={{
                          span: 5,
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
                        label="Pictures"
                        labelCol={{
                          span: 6,
                        }}>
                        <Space>
                          Left:
                          <Switch /> Right:
                          <Switch />
                        </Space>
                      </Form.Item>
                      <Form.Item
                        label="Cluster"
                        labelCol={{
                          span: 7,
                        }}>
                        {' '}
                        <Select
                          labelInValue
                          placeholder="Cluster"
                          style={{
                            width: '100%',
                          }}
                          placement="topLeft"
                          options={clusterOpsSp}
                          value={clusterCurSp}
                          onChange={(target) => {
                            setClusterCurSp({
                              value: target.value,
                              label: target.label,
                              attr: clusterOpsSp[target.value].attr,
                            })
                            toggleAnno('spConfigs')
                          }}
                        />
                      </Form.Item>
                      <Form.Item
                        label="Embedding"
                        labelCol={{
                          span: 9,
                        }}>
                        {' '}
                        <Select
                          labelInValue
                          placeholder="Embedding"
                          style={{
                            width: '100%',
                          }}
                          placement="topLeft"
                          options={embedOpsSp}
                          value={embedCurSp}
                          onChange={(value) => {
                            setEmbedCurSp(value.label)
                            toggleAnno('spConfigs')
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
                          toggleAnno('spConfigs')
                        }}>
                        Apply
                      </Button>
                    </Space>
                  </Space>
                }
                trigger="click">
                <Button type="primary" block icon={<SettingOutlined />}>
                  spConfigs
                </Button>
              </Popover>
              <Popconfirm
                title={
                  <div>
                    <div>Rename your annotation here (at most 30 characters)</div>{' '}
                    <Input
                      placeholder="Your annotations"
                      showCount
                      value={inputValue}
                      maxLength={30}
                      onChange={(event) => {
                        setInputValue(event.target.value)
                      }}
                    />
                    Previous annotation is: {nameRef.current}
                  </div>
                }
                okText="Rename"
                onConfirm={() => {
                  toggleAnno('Rename')
                }}
                cancelText="Cancel">
                <Button
                  block
                  type="primary"
                  disabled={!allowConfirm}
                  icon={<FormOutlined />}>
                  Rename
                </Button>
              </Popconfirm>
              <Popover
                title={<>Select the annotation to refine:</>}
                content={
                  <Space direction="vertical">
                    <Select
                      labelInValue
                      placeholder="Select a Refiner"
                      style={{
                        width: '100%',
                      }}
                      value={refineValue}
                      placement="topLeft"
                      options={refineOption}
                      onChange={(value) => {
                        setRefineValue(value)
                      }}
                    />
                    <Space size="small">
                      <Button
                        block
                        icon={<CloseOutlined />}
                        onClick={() => {
                          setRefineOpen(false)
                        }}>
                        Cancel
                      </Button>
                      <Button
                        block
                        icon={<ReloadOutlined />}
                        onClick={() => {
                          setRefineValue([])
                          toggleAnno('Refine_Reset')
                        }}>
                        Reset
                      </Button>
                      <Button
                        block
                        type="primary"
                        icon={<SlidersOutlined />}
                        disabled={!allowConfirm}
                        loading={loadings[0]}
                        onClick={() => {
                          enterLoading(0)
                          toggleAnno('Refine')
                        }}>
                        {loadings[0] ? "Refining" : "Refine"}
                      </Button>
                    </Space>
                  </Space>
                }
                open={refineOpen}
                onOpenChange={(newOpen) => {
                  setRefineOpen(newOpen)
                }}
                trigger="click">
                <Button block type="primary" icon={<SlidersOutlined />}>
                  Refine
                </Button>
              </Popover>
              <Popconfirm
                title={`Are you sure about \'${nameRef.current}\'?`}
                okText="Yes"
                onConfirm={() => toggleAnno('Confirm')}
                cancelText="No">
                <Button
                  block
                  type="primary"
                  disabled={!allowConfirm}
                  icon={<CheckOutlined />}>
                  Deconv
                </Button>
              </Popconfirm>
              <Popover
                title={<>Select the annotation to delete:</>}
                content={
                  <Space direction="vertical">
                    <Select
                      mode="tags"
                      placeholder="Your previous annotation"
                      style={{
                        width: '100%',
                      }}
                      value={deleteValue}
                      placement="topLeft"
                      options={nameArray}
                      onChange={(value) => {
                        setDeleteValue(value)
                      }}
                    />
                    <Space size="small">
                      <Button
                        block
                        icon={<CloseOutlined />}
                        onClick={() => {
                          setDeleteOpen(false)
                        }}>
                        Cancel
                      </Button>
                      <Button
                        block
                        icon={<ReloadOutlined />}
                        onClick={() => {
                          setDeleteValue([])
                        }}>
                        Reset
                      </Button>
                      <Button
                        block
                        danger
                        icon={<DeleteOutlined />}
                        onClick={() => {
                          toggleAnno('Delete')
                        }}>
                        Delete
                      </Button>
                    </Space>
                  </Space>
                }
                open={deleteOpen}
                onOpenChange={(newOpen) => {
                  setDeleteOpen(newOpen)
                }}
                trigger="click">
                <Button type="primary" block icon={<DeleteOutlined />}>
                  Delete
                </Button>
              </Popover>
              <Popconfirm
                title="Select the download file mode"
                okText="JSON"
                onConfirm={() => {
                  toggleAnno('Download')
                }}
                cancelText="Table">
                <Button type="primary" block icon={<CloudDownloadOutlined />}>
                  Save
                </Button>
              </Popconfirm>
            </Space>
            <div>
              <b>Current Status</b>
              <div>Mode: {brushModeState}</div>
              <div>Clustering: {clusterCurSc.label}</div>
              <div>Embedding 0: {embedCurSc}</div>
              <div>Embedding 1: {embedCurSp}</div>
              <div>Annotated Clusters: {snumRef.current}</div>
              <div>Cell number: {cellNum}</div>
            </div>
          </Space>
        </ConfigProvider>
      </Flex>
    </div>
  )
}

PairView.defaultProps = {
  height: '35rem',
  width: '55rem',
  margin: '1rem',
}

PairView.propTypes = {
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default PairView
