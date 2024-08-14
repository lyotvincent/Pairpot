import React, { useRef, useState, useEffect, useImperativeHandle } from 'react'
import * as echarts from 'echarts'
import PropTypes from 'prop-types'
import '../theme/dark'
import '../theme/vintage'
import Axis from './Axis'
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
import vega_20 from '../theme/vega_20'
import loadingTips from './LoadingTip'
import {
  ConfigProvider,
  Button,
  Space,
  Popconfirm,
  Input,
  Select,
  notification,
  Popover,
  Upload,
  Segmented,
  Slider,
  Switch,
  theme,
  Col,
  Row,
  Form,
  Flex,
  Spin,
  Tour,
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
import Loading from './Loading'
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
const {enterLoading, quitLoading} = Loading

const ScScatter = ({ scfile, spfile, setCompLoad, location, onRef, height, width, margin }) => {
  const [api, contextHolder] = notification.useNotification()
  const [isInit, setInit] = useState(false) // whether echart object is inited
  const chartRef = useRef(null) // current DOM container
  const [_data, _setData] = useState({})

  const [_scdata, setScData] = useState({})
  const [_spdata, setSpData] = useState({})
  const [_datakey, setDatakey] = useState(true)

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
  const [batchName, setBatchName] = useState("batch")
  const [batchCur, setBatchCur] = useState({})
  const [batchOps, setBatchOps] = useState([])
  const brushMode = useRef('Select')
  const [brushModeState, setBrushModeState] = useState('Select')
  const [itemSize, setItemSize] = useState(2)
  const [itemOpacity, setItemOpacity] = useState(0.8)
  const prevCluster = useRef('')
  const [clusterCur, setClusterCur] = useState({})
  const [embedCur, setEmbedCur] = useState('')
  const [embedCur1, setEmbedCur1] = useState('')
  const [clusterOps, setClusterOps] = useState({})
  const [embedOps, setEmbedOps] = useState({})
  const [xInv, setxInv] = useState(false)
  const [yInv, setyInv] = useState(false)
  const [cellNum, setCellNum] = useState(0)
  const [tourOpen, setTourOpen] = useState(false)
  const [refineOption, setRefineOption] = useState([
    { value: 2, label: 'LabelPropagation' },
  ])
  const { token } = useToken()
  const [currTip, setCurrTip] = useState(loadingTips[0])
  // ref for tours
  const UploadRef = useRef(null)
  const SwitchRef = useRef(null)
  const SettingsRef = useRef(null)
  const RenameRef = useRef(null)
  const RefineRef = useRef(null)
  const ConfirmRef = useRef(null)
  const DeleteRef = useRef(null)
  const SaveRef = useRef(null)
  const StatusRef = useRef(null)
  const steps = [
    {
      title: 'Upload File',
      description: "Click to upload your file, or drag your file here. Currently, this option accepts .h5ad file with basical attributes: .obs[\'annotation\'], and .obsm['X_umap']. These attributes can be obtained by `Scanpy` pipeline",
      target: () => UploadRef.current,
    },
    {
      title: 'Switch',
      description: 'Switch the embeddings of Single-cell data or Spatial data to display.',
      target: () => SwitchRef.current,
    },
    {
      title: 'Settings',
      description: 'Setting the appearance of LassoView, including opacity, spot-size, axis, clusters, and embeddings.',
      target: () => SettingsRef.current,
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
      title: 'Confirm',
      description: 'Confirm your activated annotation after selection and refine.',
      target: () => ConfirmRef.current,
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

  //var VisualColors = JSON.parse(JSON.stringify(vega_20)).reverse()

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
    let brushedNum = snumRef.current + 1
    let popDataset = {}
    if (doPop) {
      popDataset = __datasets.pop()
    } else {
      // default name for first brush
      nameRef.current = `annotation ${brushedNum - itemGroupRef.current.length}`
      seriesRef.current = {
        type: 'scatter',
        xAxisIndex: 1,
        yAxisIndex: 1,
        symbolSize: symbolSizeRef.current,
        name: nameRef.current,
        encode: __encode,
        large: true,
        largeThreshold: 5000,
        datasetIndex: brushedNum,
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

  const setAxis = (
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
        {
          id: 1,
          gridIndex: 1,
          name: xName1,
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
          min: Math.ceil(Math.min.apply(null, embd1_x) - 1.5),
          max: Math.ceil(Math.max.apply(null, embd1_x) + 1.5),
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
        {
          id: 1,
          gridIndex: 1,
          name: yName1,
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
          min: Math.ceil(Math.min.apply(null, embd1_y) - 1.5),
          max: Math.ceil(Math.max.apply(null, embd1_y) + 1.5),
        },
      ],
    }
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


  useImperativeHandle(onRef, () => ({  // explode trigger for parent components
    "Trigger": toggleAnno, // Trigger for useEffect
    "Tour": setTourOpen, // Open the tutorial
    "Loading": setLoadings, // set Loading status
    "Tip": setCurrTip, // set Loading Tips
  }))

  useEffect(() => {
    if (isInit) {
      // get the echart container
      var myChart = echarts.getInstanceByDom(chartRef.current)
      let _series = seriesArray
      if (commandRef.current === 'Reload') {
        setCompLoad((compLoad) => {
          let newCompLoad = { ...compLoad }
          newCompLoad['LassoView'] = true
          return newCompLoad
        })
        setCurrTip(loadingTips[1])
        enterLoading(0, setLoadings)
        SpH5adLoader(spfile).then((h5info) => {
          setSpData(h5info)
        }).catch(error => {
          console.error('Error fetching blob in LassoView:', error)
        })

        ScH5adLoader(scfile).then((h5info) => {
          console.log("All data reloaded in LassoView.")
          _setData(h5info.data)
          setTitle(h5info.title)
          setClusterOps(h5info.clusters)
          setEmbedOps(h5info.embdOps)
          setScData(h5info)
          toggleAnno("Upload")
        }).catch(error => {
          console.error('Error fetching blob in LassoView:', error)
        })
      }
      if (commandRef.current === 'Upload') {
        setCurrTip(loadingTips[2])
        let _dims = [...Object.keys(_data[0]), 'id']
        let source = _data.map((item, id) => {
          return [...Object.entries(item).map(([_, value]) => value), id]
        })
        setCellNum(source.length)
        symbolSizeRef.current = source.length > 5000 ? 2 : 4
        setItemSize(symbolSizeRef.current)
        // set annotaions or clusters
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
        let axis = setAxis(source, _dims, xName, yName, xName, yName)

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
        _series.push({
          type: 'scatter',
          symbolSize: symbolSizeRef.current,
          xAxisIndex: 1,
          yAxisIndex: 1,
          encode: {
            x: xName,
            y: yName,
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
        for (let anno of annotations) {
          // seperate the _data into subgroups by cell-types
          _snum = _snum + 1
          _series.push({
            type: 'scatter',
            symbolSize: symbolSizeRef.current,
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
        snumRef.current = _snum

        myChart.setOption(
          {
            title: [
              {
                text: title,
                left: 'center',
                textStyle: {
                  fontSize: 24,
                },
              },
            ],
            xAxis: axis.xAxis,
            yAxis: axis.yAxis,
            dataset: _datasets,
            series: _series,
            toolbox: {
              top: yInv ? '4%' : '7%',
            },
          },
          {
            replaceMerge: ['dataset', 'series'],
          }
        )
        setSeriesArray(_series)
        setBrushArray([])
        setAllowConfirm(false)
        setNameArray([])
        nameRef.current = null
        brushRef.current = []
        datasetRef.current = null
        seriesRef.current = null
        quitLoading(0, setLoadings)
      }
      if (commandRef.current === 'Setting') {
        let option = myChart.getOption()
        let _source = option.dataset[0].source
        let _dims = option.dataset[0].dimensions
        let clusterIdx = _dims.indexOf(clusterCur.label)

        // set axis
        let xName0 = embedCur + '_0'
        let yName0 = embedCur + '_1'
        let xName1 = embedCur1 + '_0'
        let yName1 = embedCur1 + '_1'
        let axis = setAxis(_source, _dims, xName0, yName0, xName1, yName1)

        // set datasets if cluster changes
        if (prevCluster.current.label !== clusterCur.label) {
          let prevAnnoLen = itemGroupRef.current.length
          let annotations = setItemGroup(_source, clusterIdx, clusterCur.attr)
          let annoLenOffset = prevAnnoLen - annotations.length
          let _datasets = setDataset(
            _source,
            _dims,
            clusterCur.label,
            annotations
          )
          let prevDatasets = option.dataset
          let myVisualMap = []
          _datasets = [..._datasets, ...prevDatasets.slice(prevAnnoLen + 1)]
          let prevSeries = option.series
          let _snum = 0
          let _series = [
            // first series are displayed in right chart
            {
              type: 'scatter',
              xAxisIndex: 1,
              yAxisIndex: 1,
              symbolSize: itemSize,
              encode: {
                x: xName1,
                y: yName1,
              },
              tooltip: { show: false }, // no tooltip
              itemStyle: {
                color: 'gray',
                opacity: itemOpacity,
              },
              large: true,
              largeThreshold: 0, //default for large mode
              datasetIndex: 0,
            },
          ]
          if (clusterCur.attr === 'categories') {
            for (let anno of annotations) {
              // seperate the _data into subgroups by cell-types
              _snum = _snum + 1
              _series.push({
                type: 'scatter',
                xAxisIndex: 0,
                yAxisIndex: 0,
                symbolSize: itemSize,
                name: anno,
                encode: {
                  // annotations are displayed in left chart
                  x: xName0,
                  y: yName0,
                  tooltip: [
                    0,
                    _dims.indexOf(xName0),
                    _dims.indexOf(yName0),
                    clusterIdx,
                    _dims[_dims.length - 1],
                  ],
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
          } else {
            _series.push({
              type: 'scatter',
              name: clusterCur.label,
              symbolSize: itemSize,
              encode: {
                // annotations are displayed in left chart
                x: xName0,
                y: yName0,
                tooltip: [
                  0,
                  _dims.indexOf(xName0),
                  _dims.indexOf(yName0),
                  clusterIdx,
                  _dims[_dims.length - 1],
                ],
              },
              itemStyle: {
                opacity: itemOpacity,
              },
              datasetIndex: _snum,
            })
            myVisualMap.push({
              show: true,
              id: clusterCur.label + '_visual',
              calculable: true,
              dimension: clusterCur.label,
              seriesIndex: 1,
              left: '10%',
              top: '5%',
              orient: 'horizontal',
              precision: 2,
              calculable: true,
              min: Math.min.apply(
                null,
                _source.map((item) => item[clusterIdx])
              ),
              max: Math.max.apply(
                null,
                _source.map((item) => item[clusterIdx])
              ),
              inRange: {
                color: ['#808080', '#FF3300'],
                opacity: [0.5, 1],
              },
              text: ['', clusterCur.label],
              textGap: 20,
              textStyle: {
                fontSize: 16,
              },
            })
          }
          // append annoated clusters
          let manualSeries = prevSeries.slice(prevAnnoLen + 1)
          for (let i = 0; i < manualSeries.length; i++) {
            manualSeries[i].datasetIndex =
              manualSeries[i].datasetIndex - annoLenOffset
            manualSeries[i].symbolSize = itemSize
            manualSeries[i].itemStyle = { opacity: itemOpacity }
            manualSeries[i].encode = {
              // manual annotations are displayed in right chart
              x: xName1,
              y: yName1,
              tooltip: [
                0,
                _dims.indexOf(xName1),
                _dims.indexOf(yName1),
                clusterIdx,
                _dims[_dims.length - 1],
              ],
            }
          }
          _series = [..._series, ...manualSeries]
          snumRef.current = _datasets.length - 1 // the sum of default and manual annotation

          myChart.setOption(
            {
              xAxis: axis.xAxis,
              yAxis: axis.yAxis,
              dataset: _datasets,
              series: _series,
              visualMap: myVisualMap,
            },
            {
              replaceMerge: ['dataset', 'series', 'visualMap'],
            }
          )
          //console.log(myChart.getOption())
          prevCluster.current = clusterCur
        } else {
          // set series only
          let _len = option.series.length
          let defAnnolen = itemGroupRef.current.length
          let _series = Array(_len)
            .fill()
            .map((_, id) => ({
              symbolSize: itemSize,
              itemStyle: { opacity: itemOpacity },
              encode: {
                x: id > 0 && id <= defAnnolen ? xName0 : xName1,
                y: id > 0 && id <= defAnnolen ? yName0 : yName1,
              },
            }))
          myChart.setOption({
            xAxis: axis.xAxis,
            yAxis: axis.yAxis,
            series: _series,
          })
          //console.log(myChart.getOption())
        }
      }
      if (commandRef.current === 'Confirm') {
        snumRef.current = snumRef.current + 1
        //console.log(`seriesNum ${snumRef.current}`)
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
        let __series = myChart.getOption().series
        __series[__series.length - 1].name = inputValue
        seriesRef.current.name = inputValue
        let prevName = nameRef.current
        nameRef.current = inputValue
        seriesRef.current = __series[__series.length - 1]
        _series.push(seriesRef.current)
        myChart.setOption({
          series: __series,
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
          .post('/api/refine', {
            data: {
              id: location.state?.st?.dataset_id,
              type: _datakey ? "sc" : "sp",
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
        jsonFile.columns = ['barcode', 'umap_0', 'umap_1', 'annotation']
        jsonFile.data = []
        for (let item of brushArray) {
          let itemData = item.data.map((index) => [
            _data[index]['index'],
            _data[index]['umap_0'],
            _data[index]['umap_1'],
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
      enterLoading(0, setLoadings)
      setCurrTip(loadingTips[0])
      // spfile.then((file) => {
      //   SpH5adLoader(file).then((h5info) => {
      //     //console.log("sp data loaded.")
      //     setSpData(h5info)
      //   })
      // }).catch(error => {
      //   console.error('Error fetching blob in LassoView:', error)
      // })
      // scfile.then((file) => {
      //   SpH5adLoader(file).then((h5info) => {
      //     console.log("All data loaded in LassoView.")
      //     _setData(h5info.data)
      //     setTitle(h5info.title)
      //     setClusterOps(h5info.clusters)
      //     setEmbedOps(h5info.embdOps)
      //     setScData(h5info)
      //     toggleAnno("Upload")
      //   })
      // }).catch(error => {
      //   console.error('Error fetching blob in LassoView:', error)
      // })
      let axis0 = Axis.setEmptyAxis(0)
      let axis1 = Axis.setEmptyAxis(1)
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
            saveAsImage: { show: true, pixelRatio: 5, },
            dataZoom: {},
          },
          iconStyle: {
            borderWidth: 1.5,
          },
          top: '7%',
        },

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
      })
      setInit(true)
      // setSeriesArray(_series)

      myChart.on('brushSelected', function (params) {
        let brushed = []
        let brushComponent = params.batch[0]
        let selectedIdx = brushComponent.selected.map((item) => {
          return item.seriesIndex
        })
        let selectedInd = brushComponent.selected.map((item) => {
          return item.dataIndex
        })
        // skip the sIdx = 0, because it's the gray scatters in index 1.
        for (let i = 0; i < selectedIdx.length; i++) {
          // search for each seriesIndex
          let sIdx = selectedIdx[i]
          let dIdx = selectedInd[i]
          for (let j = 0; j < dIdx.length; j++) {
            let dataIndex = dIdx[j]
            brushed.push(itemGroupRef.current[sIdx - 1][dataIndex])
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
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = (event) => {
        try {
          let h5info = H5adLoader(file, event, ['meta'])
          _setData(h5info.data)
          setTitle(h5info.title)
          setClusterOps(h5info.clusters)
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

  const SpH5adLoader = (file) => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.onload = (event) => {
        try {
          let h5info = H5adLoader(file, event, ['meta'])
          _setData(h5info.data)
          setTitle(h5info.title)
          setClusterOps(h5info.clusters)
          if (h5info.clusters.map((item) => item.label).includes('batch'))
            setBatchName('batch')
            let batches = [
              ...new Set(
                h5info.data.map((item) => {
                  return item['batch']
                })
              ),
            ]
            let batchIdx = Object.fromEntries(batches.map((item, idx)=>([item, idx])))
            let batch_row = Math.ceil(Math.sqrt(batches.length))
            let _row = h5info.data.map(item=>item['array_row'])
            let _col = h5info.data.map(item=>item['array_col'])
            let array_row_dist = Math.max.apply(null, _row) - Math.min.apply(null, _row)
            let array_col_dist = Math.max.apply(null, _col) - Math.min.apply(null, _col)
            h5info.data.forEach(obj=>{
              obj['array_0'] = obj['array_row']+ array_row_dist * (Math.floor(batchIdx[obj['batch']] / batch_row))
              obj['array_1'] = obj['array_col']+ array_col_dist * (batchIdx[obj['batch']] % batch_row)
            })
            h5info.embdOps = [...h5info.embdOps, { value: h5info.embdOps.length, label: "array" }]
            
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


  const JsonLoader = (file) => {
    let newTitle = file.name.replace(/\.json$/, '')
    setTitle(newTitle)
    const reader = new FileReader()
    reader.onload = (event) => {
      const jsonData = JSON.parse(event.target.result)
      _setData(jsonData)
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
      _setData(jsonData)
      toggleAnno('Upload')
    }
    reader.readAsText(file)
  }

  const beforeUpload = (file) => {
    if (file.type === 'application/json' || file.name.endsWith('.json')) {
      FileLoaderRef.current = JsonLoader
    } else if (file.name.endsWith('.h5ad')) {
      FileLoaderRef.current = ScH5adLoader
    } else {
      alert('Please select a .json or .h5ad file')
      return false
    }
    return true
  }

  const onUpload = (info) => {
    setUploading(true)
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

  return (
    <div>
      {contextHolder}
      <Flex justify="center" gap='middle'>
        <Spin spinning={loadings[0]} size="large" tip={currTip}>
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
              <div ref={UploadRef}>
                <Dragger {...upLoadProps}>
                  <p
                    className="ant-upload-drag-icon"
                    style={{ fontSize: 16, marginTop: 10, marginBottom: 10, alignItems: "center" }}>
                    <InboxOutlined />
                    <br />
                    Upload
                  </p>
                </Dragger>
              </div>
              <Segmented block ref={SwitchRef}
                options={[
                  { label: 'SC', value: true },
                  { label: 'SP', value: false },
                ]}
                value={_datakey}
                onChange={(value) => {
                  setDatakey(value)
                  let h5info = value ? _scdata : _spdata
                  _setData(h5info.data)
                  setTitle(h5info.title)
                  setClusterOps(h5info.clusters)
                  setEmbedOps(h5info.embdOps)
                  toggleAnno("Upload")
                }}
              />
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
                        span: 17,
                      }}>
                      <Form.Item label="ItemSize">
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
                        label="Inverse"
                        labelCol={{
                          span: 7,
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
                          span: 7,
                        }}>
                        <Space>
                          Left:
                          <Switch /> Right:
                          <Switch />
                        </Space>
                      </Form.Item>
                      <Form.Item
                        label="Clustering"
                        labelCol={{
                          span: 7,
                        }}>
                        {' '}
                        <Select
                          labelInValue
                          placeholder="Clustering"
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
                            toggleAnno("Settings")
                          }}
                        />
                      </Form.Item>
                      <Form.Item
                        label="Embedding 0"
                        labelCol={{
                          span: 9,
                        }}>
                        {' '}
                        <Select
                          labelInValue
                          placeholder="Embedding Left"
                          style={{
                            width: '100%',
                          }}
                          placement="topLeft"
                          options={embedOps}
                          value={embedCur}
                          onChange={(value) => {
                            setEmbedCur(value.label)
                          }}
                        />
                      </Form.Item>
                      <Form.Item
                        label="Embedding 1"
                        labelCol={{
                          span: 9,
                        }}>
                        {' '}
                        <Select
                          labelInValue
                          placeholder="Embedding Right"
                          style={{
                            width: '100%',
                          }}
                          placement="topLeft"
                          options={embedOps}
                          value={embedCur1}
                          onChange={(value) => {
                            setEmbedCur1(value.label)
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
                <Button type="primary" block icon={<SettingOutlined />} ref={SettingsRef}>
                  Settings
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
                  ref={RenameRef}
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
                          enterLoading(0, setLoadings)
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
                <Button block type="primary" icon={<SlidersOutlined />} ref={RefineRef}>
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
                  ref={ConfirmRef}
                  disabled={!allowConfirm}
                  icon={<CheckOutlined />}>
                  Confirm
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
              <b>Current Status</b>
              <div>Mode: {brushModeState}</div>
              <div>Clustering: {clusterCur.label}</div>
              <div>Embedding 0: {embedCur}</div>
              <div>Embedding 1: {embedCur1}</div>
              <div>Annotated Clusters: {snumRef.current}</div>
              <div>Cell number: {cellNum}</div>
            </div>
            <div>{JSON.stringify(_spdata[0])}</div>
          </Space>
        </ConfigProvider>
        <Tour
          open={tourOpen}
          onClose={() => setTourOpen(false)}
          steps={steps} />
      </Flex>

    </div>
  )
}

ScScatter.defaultProps = {
  height: '35rem',
  width: '55rem',
  margin: '1rem',
}

ScScatter.propTypes = {
  scfile: PropTypes.object,
  spfile: PropTypes.object,
  setCompLoad: PropTypes.func,
  location: PropTypes.object,
  onRef: PropTypes.any,
  height: PropTypes.string,
  width: PropTypes.string,
  margin: PropTypes.string,
}

export default ScScatter
