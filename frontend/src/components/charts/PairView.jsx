import React, { useRef, useState, useEffect, useImperativeHandle } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import vega_20 from '../theme/vega_20'
import strokeColor from '../theme/strokeColor'
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
import Dataset from './Dataset'
import Loading from './Loading'
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
  Spin,
  Tour,
  Progress
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
  ContactsOutlined,
} from '@ant-design/icons'
import loadingTips from './LoadingTip'
echarts.use([
  GraphicComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  LineChart,
  CanvasRenderer,
  UniversalTransition,
])

const { useToken } = theme
const { quitLoading, enterLoading } = Loading

const PairView = ({ spfile, scfile, setCompLoad, onRef, height, width, margin, meta, progress }) => {
  const [api, contextHolder] = notification.useNotification()
  const [isInit, setInit] = useState(false) // whether echart object is inited
  const chartRef = useRef(null) // current DOM container
  const [_scdata, _setScData] = useState({})
  const [_spdata, _setSpData] = useState({})
  const [seriesArray, setSeriesArray] = useState([])
  const [nameArray, setNameArray] = useState([])
  const snumRef = useRef()     // current series nums
  const datasetRef = useRef()  // current dataset
  const seriesRef = useRef()   // current series
  const commandRef = useRef()  // current command
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
  const [title, setTitle] = useState('PairView Title')
  const symbolSizeRef = useRef('')
  const [loadings, setLoadings] = useState([false, false])
  const brushMode = useRef('Select')
  const [brushModeState, setBrushModeState] = useState('Select')
  const [currProps, setCurrProps] = useState([])

  // Item size and Opacity for Single-cell and Spatial
  const [itemSizeSc, setItemSizeSc] = useState(2)
  const [itemOpacitySc, setItemOpacitySc] = useState(0.8)
  const [itemSizeSp, setItemSizeSp] = useState(4)
  const [itemOpacitySp, setItemOpacitySp] = useState(0.8)
  const prevCluster = useRef('')

  // Clusters and Embeddings for Single-cell
  const [clusterCurSc, setClusterCurSc] = useState({})
  const [embedCurSc, setEmbedCurSc] = useState('')
  const [clusterOpsSc, setClusterOpsSc] = useState({})
  const [embedOpsSc, setEmbedOpsSc] = useState({})
  const [annoCurSc, setAnnoCurSc] = useState([])

  // Clusters, Batches and Embeddings for SRT
  const [clusterCurSp, setClusterCurSp] = useState({})
  const [embedCurSp, setEmbedCurSp] = useState('')
  const [clusterOpsSp, setClusterOpsSp] = useState([])
  const [embedOpsSp, setEmbedOpsSp] = useState([])
  const [annoCurSp, setAnnoCurSp] = useState([])
  const stAnnoLength = useRef(0)
  const [batchName, setBatchName] = useState("batch")
  const [batchCur, setBatchCur] = useState({})
  const [batchOps, setBatchOps] = useState([])
  const [propOpsSp, setPropOpsSp] = useState([])
  const [propCurSp, setPropCurSp] = useState([])

  const [tourOpen, setTourOpen] = useState(false)
  const [currTip, setCurrTip] = useState(loadingTips[0])

  // ref for tours
  const scConfigsRef = useRef(null)
  const spConfigsRef = useRef(null)
  const RenameRef = useRef(null)
  const RefineRef = useRef(null)
  const DeconvRef = useRef(null)
  const DeleteRef = useRef(null)
  const SaveRef = useRef(null)
  const StatusRef = useRef(null)
  const steps = [
    {
      title: 'scConfigs',
      description: 'Setting the appearance of LassoView.',
      target: () => scConfigsRef.current,
    },
    {
      title: 'spConfigs',
      description: 'Setting the appearance of LassoView.',
      target: () => spConfigsRef.current,
    },
    {
      title: 'Rename',
      description: 'Rename the activated lasso regions. The lasso button is at the topRight of the Charts. You can lasso spots or cells in the left chart, and the results would display in the right chart.',
      target: () => RenameRef.current,
    },
    {
      title: 'Refine',
      description: 'Select a refiner to Refine the acivated lasso regions. The refined results would display in the right chart.',
      target: () => RefineRef.current,
    },
    {
      title: 'Deconv',
      description: 'Confirm your activated annotation after selection and refine.',
      target: () => DeconvRef.current,
    },
    {
      title: 'Delete',
      description: 'Delete your confirmed annotation.',
      target: () => DeleteRef.current,
    },
    {
      title: 'Save',
      description: 'Save all your online-defined annotations as a json file.',
      target: () => SaveRef.current,
    },
    {
      title: 'Current Status',
      description: 'The Status of Current LassoView, including current Lasso mode, clustering, embeddings, cell numbers and annoatated clusters.',
      target: () => StatusRef.current,
    },
  ]


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
        nameLocation: 'center',
        nameGap: 30,
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
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = (event) => {
        try {
          let h5info = H5adLoader(file, event, ['meta'])
          _setSpData(h5info.data)
          setTitle(meta?.st.dataset_id)
          if (h5info.clusters.map((item) => item.label).includes('batch'))
            setBatchName('batch')
          setClusterOpsSp(h5info.clusters)
          setEmbedOpsSp(h5info.embdOps)
          setPropOpsSp(h5info.propOps)
          return resolve({
            data: h5info.data,
            title: h5info.title,
            clusters: h5info.clusters,
            embdOps: h5info.embdOps,
            propOps: h5info.propOps,
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

  const ScH5adLoader = (file) => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = (event) => {
        try {
          let h5info = H5adLoader(file, event, ['meta'])
          _setScData(h5info.data)
          setClusterOpsSc(h5info.clusters)
          setEmbedOpsSc(h5info.embdOps)
          resolve({
            data: h5info.data,
            clusters: h5info.clusters,
            embdOps: h5info.embdOps
          })
        } catch (error) {
          reject(error)
        }
      }
      reader.onloadend = () => {
        //console.log("Load sc data finished.")
      }
      reader.onerror = (error) => {
        reject(error)
      }
      reader.readAsArrayBuffer(file)
    })

  }
  const actionRef = useRef(0)
  const toggleAnno = (command) => {
    commandRef.current = command
    actionRef.current = actionRef.current + 1
    setAction(actionRef.current)
  }

  useImperativeHandle(onRef, () => ({  // explode trigger for parent components
    "Trigger": toggleAnno, // Trigger for useEffect
    "Tour": setTourOpen, // Open the tutorial
    "Loading": setLoadings, // set Loading status
    "Tip": setCurrTip, // set Loading Tips
  }), [])

  useEffect(() => {
    if (isInit) {
      // get the echart container
      var myChart = echarts.getInstanceByDom(chartRef.current)
      let _series = seriesArray
      if (commandRef.current === 'Reload') {
        setCompLoad((compLoad) => {
          let newCompLoad = { ...compLoad }
          newCompLoad['PairView'] = true
          return newCompLoad
        })
        enterLoading(1, setLoadings)
        setCurrTip(loadingTips[1])
        Promise.all([SpH5adLoader(spfile), ScH5adLoader(scfile)]).then(() => {
          console.log("All data reloaded in PairView.")
          toggleAnno("Upload")
        }).catch(error => {
          console.error('Error fetching blob in PairView:', error)
        })
      }
      if (commandRef.current === 'Upload') {
        enterLoading(1, setLoadings)
        setCurrTip(loadingTips[2])
        let _scdims = [...Object.keys(_scdata[0]), 'id']
        let _scsource = _scdata.map((item, id) => {  // 2d array
          return [...Object.entries(item).map(([_, value]) => value), id]
        })

        let _spdims = [...Object.keys(_spdata[0]), 'id']
        let _spsource = _spdata.map((item, id) => {  // 2d array
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
        //setClusterOpsSc([{ value: 0, label: scdefaultAnno, attr: 'categories' }, { value: 1, label: 'umap_0', attr: 'numeric' }])
        let _scannotations = setItemGroup(_scsource, _scdims.indexOf(scdefaultAnno), 'categories', true)
        setAnnoCurSc(_scannotations)

        // set sp annotations
        let spdefaultAnno
        if (_spdims.includes('annotation')) {
          spdefaultAnno = 'annotation'
        } else {
          spdefaultAnno = 'leiden-1'
        }
        setClusterCurSp({ value: 0, label: spdefaultAnno, attr: 'categories' })
        let _spannotations = setItemGroup(_spsource, _spdims.indexOf(spdefaultAnno), 'categories', false)
        setAnnoCurSp(_spannotations)
        stAnnoLength.current = _spannotations.length
        setPropCurSp({ value: 0, label: "Cell2Location" })

        // set sp batches
        let batches = []
        if (_spdims.includes(batchName)) {
          batches = setItemGroup(_spsource, _spdims.indexOf(batchName), 'categories', false)
          setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
          setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))
        }
        else {
          setBatchName(null)
          batches = ['batch 1']
          setBatchCur({ value: 0, label: batches[0], attr: 'categories' })
          setBatchOps(batches.map((item, id) => ({ value: id, label: item, attr: 'categories' })))
        }

        // set sc embeddings
        let _scaxis = setAxis(_scsource, _scdims, 'X_umap_0', 'X_umap_1', 0)
        setEmbedCurSc('X_umap')

        // set sp embeddings
        let _spaxis = setAxis(_spsource, _spdims, 'spatial_0', 'spatial_1', 1)
        setEmbedCurSp('spatial')


        // set datasets
        let _scdatasets = Dataset.setDataset(_scsource, _scdims, scdefaultAnno, _scannotations, 0)
        let _spdatasets = Dataset.setBatchDataset(_spsource, _spdims, spdefaultAnno, _spannotations, batches[0], batchName, _scdatasets.length, false)
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
              x: 'X_umap_0',
              y: 'X_umap_1',
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
              x: 'spatial_0',
              y: 'spatial_1',
              itemName: anno,
            },
            emphasis: {
              focus: 'series',
            },
            datasetIndex: _snum,
          })
        }
        snumRef.current = _snum
        brushRef.current = []
        // set visualMap
        myChart.setOption({
          color: vega_20,
          title: [
            {
              text: "Single-cell Data",
              left: '18%',
              top: '10%',
              textStyle: {
                fontSize: 20,
                fontWeight: 'normal',
              },
            },
            {
              text: "Spatial Omics Data",
              right: '9%',
              top: '10%',
              textStyle: {
                fontSize: 20,
                fontWeight: 'normal',
              },
            },
            {
              text: title,
              left: 'center',
              textStyle: {
                fontSize: 24,
              },
            },
          ],
          toolbox: {
            feature: {
              saveAsImage: {
                show: true,
                pixelRatio: 4,
                name: `PairView_${title}`,
              },
            }
          },
          xAxis: [_scaxis.xAxis, _spaxis.xAxis],
          yAxis: [_scaxis.yAxis, _spaxis.yAxis],
          grid: [
            {
              top: '18%',
              left: '8%',
              width: '42%',
              bottom: '18%',
            },
            {
              top: '18%',
              width: '42%',
              right: '1%',
              bottom: '18%',
            },
          ],
          dataset: _datasets,
          series: _series,
          legend: {
            type: 'scroll',
            animation: false,
            pageIconColor: token.colorPrimaryActive,
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

        myChart.on('click', function (params) {
          if (params.componentType === 'series' &&
            params.componentSubType === 'scatter'
          ) {
            hoverClu(params.seriesName)
          }
        })
        quitLoading(1, setLoadings)
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
        if (stAnnoLength.current === 0) {
          prevAnnoLen = prevAnnoLen + 1
        }
        let annoLenOffset = _sclen - 1 - annotations.length
        setAnnoCurSc(annotations)
        if (annotations.length === 0) {
          annoLenOffset = annoLenOffset - 1
        }
        let _scdatasets = Dataset.setDataset(_scsource, _scdims, clu, annotations, 0)
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
        //console.log(myChart.getOption())
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
        if (embedCurSp === 'spatial' && _spdims.includes("array_row") && _spdims.includes("array_col")) {
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
        let attr = clusterCurSp.attr
        // use cell proportions from current dcv method.
        if (annoCurSc.includes(clu)) {
          let _clu = clusterOpsSp.find(item => item.label === `${propCurSp.label}_${annoCurSc.indexOf(clu)}`)
          setClusterCurSp(_clu)
          clu = _clu.label
          attr = _clu.attr
        }
        let cluIdx = _spdims.indexOf(clu)
        let prevAnnoLen = stAnnoLength.current + _sclen + 1
        if (stAnnoLength.current === 0) {
          prevAnnoLen = prevAnnoLen + 1
        }
        let annotations = setItemGroup(_spsource, cluIdx, attr, false)  // spconfig do not change the item group
        let annoLenOffset = stAnnoLength.current - annotations.length
        if (stAnnoLength.current === 0) {
          annoLenOffset = annoLenOffset + 1
        }
        stAnnoLength.current = annotations.length
        setAnnoCurSp(annotations)
        let _spdatasets = Dataset.setBatchDataset(_spsource, _spdims, clu, annotations, batchCur.label, batchName, _sclen, false)
        let _datasets = [...option.dataset.slice(0, _sclen), ..._spdatasets, ...option.dataset.slice(prevAnnoLen, prevAnnoLen + 1)]
        //console.log(prevAnnoLen, _datasets)
        // set series and VisualMap
        let _spSeries = []
        let _snum = _sclen
        let myVisualMap = []
        if (attr === 'categories') {
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
            datasetIndex: _snum + 1,
          })

          myVisualMap.push({
            show: true,
            id: clu + '_sp',
            calculable: true,
            dimension: clu,
            seriesIndex: _sclen - 1,
            left: '60%',
            top: '0%',
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
        //console.log(manualSeries)
        let _series = [...option.series.slice(0, _sclen - 1), ..._spSeries, ...manualSeries]
        snumRef.current = _series.length - manualSeries.length
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
        //console.log("sp", _visualMap)

        myChart.setOption({
          xAxis: [option.xAxis[0], spaxis.xAxis],
          yAxis: [option.yAxis[0], spaxis.yAxis],
          series: _series,
          dataset: _datasets,
          visualMap: _visualMap,
        },
          {
            replaceMerge: ['series', 'dataset', 'visualMap'],
          })
        //console.log(myChart.getOption())
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
      if (commandRef.current === 'Deconv') {
        let starttime = Date.now()
        axios
          .post('/api/deconv', {
            data: {
              id: meta.st?.dataset_id,
              anno: brushRef.current,
              starttime: starttime,
            },
            headers: {
              'Content-Type': 'application/json',
            },
          }).then((response) => {
            if (response.data.success) {
              setBrushArray([  // set current brush array
                ...brushArray,
                {
                  name: nameRef.current,
                  data: brushRef.current,
                  value: snumRef.current - itemGroupRef.current.length,
                },
              ])
              let endtime = response.data.endtime
              let props = response.data.props
              let option = myChart.getOption()   // set options
              // set sp source
              let _sclen = annoCurSc.length + 1
              if (annoCurSc.length === 0) {
                _sclen = _sclen + 1
              }
              let _dataset = option.dataset
              let _spsource = _dataset[_sclen].source
              let _spdims = _dataset[_sclen].dimensions
              // set current customized props
              setCurrProps(_spsource.map((item, id) => [item[_spdims.indexOf("index")], props[id]]))
              setCurrProps(_spsource.map((item, id) => [item[_spdims.indexOf("index")], props[id]]))
              if (_spdims[_spdims.length - 1] === "Customized cell props") { // if exists Cell props
                _spsource = _spsource.map((item, id) => {  // 2d array
                  item[item.length - 1] = props[id]
                  return item
                })
                let _clusterOpsSp = clusterOpsSp
                _clusterOpsSp[clusterOpsSp.length - 1] = { value: clusterOpsSp.length - 1, label: "Customized cell props", attr: 'values' }
                setClusterOpsSp(_clusterOpsSp)
                setClusterCurSp({ value: clusterOpsSp.length - 1, label: "Customized cell props", attr: 'values' })
              } else {
                _spdims = [..._spdims, "Customized cell props"]
                _spsource = _spsource.map((item, id) => {  // 2d array
                  return [...item, props[id]]
                })
                setClusterOpsSp([...clusterOpsSp, { value: clusterOpsSp.length, label: "Customized cell props", attr: 'values' }])
                setClusterCurSp({ value: clusterOpsSp.length, label: "Customized cell props", attr: 'values' })
              }

              // set Options
              _dataset[_sclen] = {
                dimensions: _spdims,
                source: _spsource,
              }
              myChart.setOption({
                dataset: _dataset,
              })
              api.success({
                message: `Annotation Refined in ${(endtime - starttime) / 1000}s`,
                description: response.data.message,
                placement: 'topRight',
              })
              toggleAnno("spConfigs")
            } else {
              let endtime = response.data.endtime
              api.warning({
                message: `Annotation does not Refined in ${(endtime - starttime) / 1000}s`,
                description: response.data.message,
                placement: 'topRight',
              })
            }
          }).catch((error) => {
            api.error({
              message: `Annotation failed for ${error}`,
              placement: 'topRight',
            })
          })
      }
      if (commandRef.current === 'Rename') {
        let option = myChart.getOption()
        let _series = option.series
        _series[_series.length - 1].name = inputValue
        let prevName = nameRef.current
        nameRef.current = inputValue
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
        // modify corresponding item in brush array
        let _newbrush = brushArray
        let _idx = _newbrush.findIndex(item => item.name === prevName)
        _newbrush[_idx].name = inputValue
        setBrushArray(_newbrush)
      }
      if (commandRef.current === 'Refine') {
        let starttime = Date.now()
        axios
          .post('/api/refine', {
            data: {
              id: meta.st?.dataset_id,
              type: "sc",
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
              let option = myChart.getOption()  // set options
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
        quitLoading(0, setLoadings)
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
        //console.log(myChart.getOption())
        api.warning({
          message: 'Annotation Deleted',
          description:
            'Previous annotation is deleted, but you can recover them until leaving this page.',
          placement: 'topRight',
        })
      }
      if (commandRef.current === 'Download') {
        let jsonFile = {}
        jsonFile.sc_columns = ['barcode', 'annotation']
        jsonFile.sc_annotation = []
        for (let item of brushArray) {
          let itemData = item.data.map((index) => [
            _scdata[index]['index'],
            item.name,
          ])
          jsonFile.sc_annotation = [...jsonFile.sc_annotation, ...itemData]
        }
        jsonFile.sp_columns = ['barcode', 'proportion']
        jsonFile.sp_proportion = currProps
        const fileStream = JSON.stringify(jsonFile)
        const blob = new Blob([fileStream])
        saveAs(blob, `Annotation-${new Date().getTime()}.json`)
      }
    } else {
      // init the echart container
      var myChart = echarts.init(chartRef.current)
      enterLoading(1, setLoadings)
      setCurrTip(loadingTips[0])
      // Promise.all([spfile.then((file) => SpH5adLoader(file)),
      // scfile.then((file) => ScH5adLoader(file))]).then(() => {
      //   console.log("All data loaded in PairView.")
      //   toggleAnno("Upload")
      // }).catch(error => {
      //   console.error('Error fetching blob in PairView:', error)
      // })

      let _scaxis = setEmptyAxis(0)
      let _spaxis = setEmptyAxis(1)
      myChart.setOption({
        textStyle: {
          fontFamily: 'Arial',
        },
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
            saveAsImage: {
              show: true,
              pixelRatio: 3,
              name: `PairView_${title}`,

            },
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
            left: '8%',
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
        let selected = brushComponent.selected
        console.log(brushRef.current)
        if (brushRef.current.length > 0) {  // existing brushed series
          console.log(selected.length, snumRef.current)
          if (selected.length > snumRef.current) {
            selected = selected.slice(0, -1)  // brushed series is rendered
          }
        }
        let selectedIdx = selected.map((item) => {
          return item.seriesIndex
        })
        let selectedInd = selected.map((item) => {
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

  const hoverClu = (label) => {
    let _clu = clusterOpsSp.find((item) => item.label === label)
    if (typeof _clu !== 'undefined') {
      setClusterCurSp(_clu)
      toggleAnno("spConfigs")
    }
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
        <Spin
          spinning={loadings[1]}
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
                      <Button icon={<ReloadOutlined />}
                        onClick={() => {
                          setItemSizeSc(4)
                          setItemOpacitySc(0.8)
                          setBrushModeState('Select')
                          setyInv(false)
                          setxInv(false)
                          toggleAnno("scConfigs")
                        }}
                      >Reset</Button>
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
                <Button type="primary" block icon={<SettingOutlined />} ref={scConfigsRef}>
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
                        label="Props"
                        labelCol={{
                          span: 6,
                        }}>
                        <Select
                          labelInValue
                          placeholder="Methods"
                          style={{
                            width: '100%',
                          }}
                          placement="topLeft"
                          options={propOpsSp}
                          value={propCurSp}
                          onChange={(target) => {
                            setPropCurSp({
                              value: target.value,
                              label: target.label,
                            })
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
                        label="Reduction"
                        labelCol={{
                          span: 7,
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
                            toggleAnno("spConfigs")
                          }}
                        />
                      </Form.Item>
                    </Form>
                    <Space size="small">
                      <Button
                        icon={<ReloadOutlined />}
                        onClick={() => {
                          setItemSizeSp(4)
                          setItemOpacitySp(0.8)
                          setBrushModeState('Select')
                          setyInv(false)
                          setxInv(false)
                          toggleAnno("spConfigs")
                        }}
                      >Reset</Button>
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
                <Button type="primary" block icon={<SettingOutlined />} ref={spConfigsRef}>
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
                  icon={<FormOutlined />}
                  ref={RenameRef}>
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
                          enterLoading(0, setLoadings)
                          toggleAnno('Refine')
                        }}
                      >
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
                <Button block type="primary" icon={<SlidersOutlined />} ref={RefineRef}>
                  Refine
                </Button>
              </Popover>
              <Popconfirm
                title={`Are you sure about \'${nameRef.current}\'?`}
                okText="Yes"
                onConfirm={() => toggleAnno('Deconv')}
                cancelText="No">
                <Button
                  block
                  type="primary"
                  disabled={!allowConfirm}
                  icon={<CheckOutlined />}
                  ref={DeconvRef}>
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
                <Button type="primary" block icon={<DeleteOutlined />} ref={DeleteRef}>
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
                <Button type="primary" block icon={<CloudDownloadOutlined />} ref={SaveRef}>
                  Save
                </Button>
              </Popconfirm>
            </Space>
            <div ref={StatusRef}>
              <div>
                <b>Current Status</b>
                <div>Mode: {brushModeState}</div>
                <div>Clustering: {clusterCurSc.label}</div>
                <div>Embedding 0: {embedCurSc}</div>
                <div>Embedding 1: {embedCurSp}</div>
                <div>Deconv Method: {propCurSp.label}</div>
                <div>Annotated Clusters: {snumRef.current}</div>
                <div>Cell number: {cellNum}</div>
              </div>
            </div>
            {/* <div>{JSON.stringify(brushArray)}</div> */}
          </Space>
        </ConfigProvider>
        <Tour open={tourOpen} onClose={() => setTourOpen(false)} steps={steps} />
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
  spfile: PropTypes.object,
  scfile: PropTypes.object,
  setCompLoad: PropTypes.func,
  onRef: PropTypes.any,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
  meta: PropTypes.any,
  progress: PropTypes.number,
}

export default PairView
