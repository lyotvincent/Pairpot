import React, { useState } from 'react'
import { UserOutlined, DownOutlined } from '@ant-design/icons'
import { AutoComplete, Input, Col, Row, Popover, ConfigProvider, Button } from 'antd'
import { useNavigate, Link } from 'react-router-dom'
import Loading from '../charts/Loading'
import axios from 'axios'
import TextCollapse from './TextCollapse'
import { useMutation } from 'react-query'
const { enterLoading, quitLoading } = Loading
const renderTitle = (title) => (
  <span>
    {title}
    <Link
      style={{
        float: 'right',
      }}
      to="/database"
    >
      more
    </Link>
  </span>
)
const renderItem = (title, key) => ({
  value: key,
  label: (
    <div
      style={{
        display: 'flex',
        justifyContent: 'space-between',
      }}
    >
      <Popover content={title} zIndex={2000}>
        {title.length > 90 ? `${title.substring(0, 90)}...` : title}
      </Popover>
      <span>
        <UserOutlined /> {key}
      </span>
    </div>
  ),
})
const options = [
  {
    label: renderTitle('Spatial Transcriptomics Studies'),
    options: [renderItem('Spatial transcriptomics map of the embryonic mouse brain – a tool to explore neurogenesis', 'STDS0000235'),
    renderItem("SAW: An efficient and accurate data analysis workflow for Stereo-seq spatial transcriptomics]", 'STDS0000234'),
    renderItem("Integration of spatial and single-cell transcriptomics localizes epithelial cell-immune cross-talk in kidney injury", 'STDS0000007'),
    renderItem('Spatial transcriptomics combined with single-cell RNA-sequencing unravels the complex inflammatory cell network in atopic dermatitis [ST]', 'STDS0000212'),
    renderItem("Spatial transcriptomics of de novo NEPC and ARPC", 'STDS0000227'),
    renderItem("Tertiary lymphoid structures generate and propagate anti-tumor antibody-producing plasma cells in renal cell cancer", 'STDS0000223'),
    renderItem("Spatial resolution of cellular senescence dynamics in colorectal liver metastasis", 'STDS0000219'),
    renderItem("Spatial transcriptomics of adenoid cystic carcinoma of the lacrimal gland", 'STDS0000221'),
    renderItem("Spatial transcriptomics (Visium, 10x Genomics) data of Duchenne mouse models", 'STDS0000204'),
    renderItem("Discovering Haematoma-Stimulated Circuits for Secondary Brain Injury after Intraventricular Haemorrhage by Spatial Transcriptome Analysis", 'STDS0000201'),
    renderItem("Moxibustion Improves Hypothalamus Aqp4 Polarization in APP/PS1 Mice: Evidence from Spatial Transcriptomics", 'STDS0000189'),
    renderItem("A topographic atlas defines developmental origins of cell heterogeneity in the human embryonic lung [ST]", "STDS0000175"),
    renderItem("Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment", 'STDS0000153'),
    renderItem("Intratumor heterogeneity and T cell exhaustion in primary CNS lymphoma", 'STDS0000152'),
    renderItem("Spatially resolved transcriptomics revised human distal lung epithelial hierarchy [spatial transcriptomics]", 'STDS0000114'),
    renderItem("Spatial detection of fetal marker genes expressed at low level in adult human heart tissue", 'STDS0000014'),
    renderItem("A spatially resolved brain region- and cell type-specific isoform atlas of the postnatal mouse brain", "STDS0000004"),
    renderItem("A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart", 'STDS0000015'),
    renderItem("Mouse Brain Serial Section 1 (Sagittal-Anterior)", 'STDS0000018'),
    renderItem("Mouse Brain Serial Section 2 (Sagittal-Anterior)", 'STDS0000019'),
    renderItem("Mouse Brain Serial Section 2 (Sagittal-Posterior)", 'STDS0000020'),
    renderItem("Mouse Brain Serial Section 1 (Sagittal-Posterior)", 'STDS0000021'),
    renderItem("Mouse Brain Section (Coronal)", 'STDS0000022'),
    renderItem("Mouse Kidney Section (Coronal)", 'STDS0000023'),
    renderItem("Human Lymph Node", 'STDS0000024'),
    renderItem("Human Heart", 'STDS0000025'),
    renderItem("Human Breast Cancer (Block A Section 2)", 'STDS0000026'),
    renderItem("Human Breast Cancer (Block A Section 1)", 'STDS0000027'),
    renderItem("Invasive Ductal Carcinoma Stained With Fluorescent CD3 Antibody", 'STDS0000028'),
    renderItem("Adult Mouse Brain Section 1 (Coronal). Stains: DAPI, Anti-NeuN", 'STDS0000029'),
    renderItem("Adult Mouse Brain Section 2 (Coronal). Stains: DAPI, Anti-GFAP, Anti-NeuN", 'STDS0000030'),

    renderItem("Human Breast Cancer: Whole Transcriptome Analysis", 'STDS0000038'),
    renderItem("Human Breast Cancer: Targeted, Immunology Panel", 'STDS0000039'),
    renderItem("Human Prostate Cancer, Adenocarcinoma with Invasive Carcinoma (FFPE)", 'STDS0000046'),
    renderItem("Human Prostate Cancer, Adjacent Normal Section with IF staining (FFPE)", 'STDS0000047'),
    renderItem("Normal Human Prostate (FFPE)", 'STDS0000048'),
    renderItem("Human Breast Cancer: Ductal Carcinoma In Situ, Invasive Carcinoma (FFPE)", 'STDS0000049'),
    renderItem("Adult Mouse Brain, IF stained (FFPE)", 'STDS0000050'),

    renderItem("Adult Mouse Kidney (FFPE)", 'STDS0000051'),
    renderItem("Adult Mouse Brain (FFPE)", 'STDS0000052'),
    renderItem("Spatial maps of prostate cancer transcriptomes reveal an unexplored landscape of heterogeneity", 'STDS0000053'),
    renderItem("10x Visium Spatial transcriptomics of murine colon at d14 (mucosa healing) in B cell sufficient/deficient mice", 'STDS0000063'),
    renderItem("A new epithelial cell subpopulation predicts response to surgery, chemotherapy, and immunotherapy in bladder cancer", 'STDS0000067'),
    renderItem("Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas 9606", 'STDS0000073'),
    renderItem("Cell types and clonal relations in the mouse brain revealed by single-cell and spatial transcriptomics", 'STDS0000085'),
    renderItem("Interleukin-17 governs hypoxic adaptation of injured epithelium [spatial transcriptomics]", 'STDS0000089'),
    renderItem("Spatial transcriptome analysis of gastric cancer which GAN-KP transplanted C57BL/6J mouse", 'STDS0000093'),

    renderItem("The spatial landscape of gene expression isoforms in tissue sections", 'STDS0000107'),
    renderItem("Human Cervical Cancer (FFPE)", 'STDS0000113'),

    renderItem("Dietary palmitic acid promotes a prometastatic memory based on Schwann Cell activation", 'STDS0000120'),
    renderItem("Spatially resolved transcriptomic analysis of acute kidney injury in a female murine model", 'STDS0000121'),
    renderItem("Characterizing Neonatal Heart Maturation, Regeneration, and Scar Resolution Using Spatial Transcriptomics [spatial]", 'STDS0000124'),
    renderItem("Human Prostate Cancer, Acinar Cell Carcinoma (FFPE)", 'STDS0000131'),

    renderItem("Spatial transcriptomics: The effect of consecutive slices data integration on accurate cell type annotation and clustering", 'STDS0000224')

    ],
  },
]

const Search = ({ setLocState, setCompLoad, setSelectedKey }) => {
  const [loadText, setLoadText] = useState("Browse")
  const navigate = useNavigate()
  const [loading, setLoading] = useState([])
  const [open, setOpen] = useState(false)
  const [status, setStatus] = useState("Success")
  const [searchText, setSearchText] = useState("")
  const [MyOptions, setOptions] = useState(options)


  // 下拉框根据输入来动态更新
  const handleSearch = (value) => {
    setSearchText(value) // 更新搜索框状态
    if (value.length > 0) {
      const filteredOptions = options.map(optionGroup => ({
        ...optionGroup,
        options: optionGroup.options.filter(item => item.label.props.children[0].props.children.toLowerCase().includes(value.toLowerCase()))
      }))
      setOptions(filteredOptions)
    } else {
      setOptions(options)  // 如果没有输入内容，则重置为初始的选项
    }
  }

  const fetch = useMutation({
    mutationKey: ['example'],
    mutationFn: (e) => {
      axios({
        method: 'GET',
        url: '/api/example',
        params: {
          id: e
        },
      }).then((response) => {
        let dataCol = response.data.attributes
        let spitem = response.data.data[0]
        let values = Object.fromEntries(
          dataCol.map((k, i) => [k, spitem[i]])
        )
        let scitem = response.data.data[1]
        let state = {
          st: values
        }
        if (typeof scitem !== 'undefined') {
          let scvalues = Object.fromEntries(
            dataCol.map((k, i) => [k, scitem[i]])
          )
          state['sc'] = scvalues
        }
        quitLoading(0, setLoading)
        setLoadText("Browse")
        setLocState(state)
        setCompLoad({  // reset ComponentLoad
          "LassoView": false,
          "LayerView": false,
          "PairView": false,
          "MarkerTable": false,
          "NetWorkRelation": false,
          "CPDBHeatmap": false,
        })
        setSelectedKey("MetaInfo")
      })
    }
  })
  return (
    <ConfigProvider theme={{
      components: {
        Popover: {
          colorBgElevated: '#ddf7d4',
        }
      }
    }}>
      <Row gutter={8}>
        <Col flex="auto">
          <AutoComplete
            popupClassName="certain-category-search-dropdown"
            popupMatchSelectWidth={800}
            style={{
              width: "100%"
            }}
            options={MyOptions}
            open={open}
            onBlur={() => { setOpen(false) }}
            onSelect={(value) => {
              setOpen(false)
              setSearchText(value)
              // fetch
              let e = value
              setSearchText(e)
              setOpen(false)
              setStatus("Success")
              enterLoading(0, setLoading)
              setLoadText("Fetching data...")
              fetch.mutate(e)
            }}
            onSearch={handleSearch}
            size="large"
          >
            {/* <Input.Search */}
            <Input
              size="large"
              placeholder="Select a Study to Browse."
              // enterButton={loadText}
              status={status}
              loading={loading[0]}
              onClick={() => { setOpen(!open) }}
              onChange={(e) => setSearchText(e.target.value)} // input on time
              suffix={<DownOutlined style={{ cursor: 'pointer' }} onClick={() => { setOpen(!open) }} />}
            // onSearch={(e) => {
            //   // setSearchText(e)
            //   if (e.length > 0) {
            //     if(!e.startsWith("STDS00")){
            //       // setLoadText("Browse")
            //       // setStatus("error")
            //       let text = MyOptions[0].options.length == 0 ? "STDS0000235" : MyOptions[0].options[0].value
            //       setSearchText(text)
            //       setOpen(false)
            //       setStatus("Success")
            //       enterLoading(0, setLoading)
            //       setLoadText("Fetching data...")
            //       fetch.mutate(text)
            //     }
            //     else{
            //       setSearchText(e)
            //       setOpen(false)
            //       setStatus("Success")
            //       enterLoading(0, setLoading)
            //       setLoadText("Fetching data...")
            //       fetch.mutate(e)
            //     }
            //   } else {
            //     setStatus("error")
            //   }
            // }} 
            />
          </AutoComplete>
        </Col>

        <Col>
          <Button
            onClick={() => {
              // const myKey = searchText
              navigate('/database', { state: { searchText } })
            }}
          >
            {"Browse in Database >>>"}
          </Button>
        </Col>
      </Row>
    </ConfigProvider>
  )
}

export default Search
