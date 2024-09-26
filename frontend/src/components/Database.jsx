import React, { useState, useEffect, View } from 'react'
// import TrackVisibility from 'react-on-screen'
import 'bootstrap/dist/css/bootstrap.min.css'
// import { ReactComponent as Female } from '../assets/img/female.svg'
// import { ReactComponent as FemaleOrgans } from '../assets/img/female-organs.svg'
//import SideMenu from './utils/SideMenu'
import Species from './metas/Species'
import Technologies from './metas/Technologies'
import Tissues from './metas/Tissues'
import Topics from './metas/Topics'
import axios from 'axios'
import { useQuery } from 'react-query'
import { Card, Layout, theme, Tag, Row, Col, Space, Divider, ConfigProvider } from 'antd'
import logoFig from "../assets/img/mylogo.png"
import TagCollapse from './utils/TagCollapse'
// import loadingTips from './charts/LoadingTip'
import Loading from './charts/Loading'
import contactImg from "../assets/img/contact.png"
import {
  DatabaseOutlined,
  ProfileOutlined,
  OneToOneOutlined,
} from '@ant-design/icons'
import DatasetTab from './utils/DatasetTab'
import Search from './utils/SearchEngine'
import SampleTab from './utils/SampleTab'
const { Content, Footer } = Layout
const { enterLoading } = Loading

export const Database = () => {
  //const [collapsed, setCollapsed] = useState(false)
  const [srcMeta, setSrcMeta] = useState({}) // get the meta data from childrens
  const [src, setSrc] = useState([1, 2, 3])
  const [selectedTags, setSelectedTags] = React.useState(null)
  const [searchValue, steSearchValue] = useState('') // 搜索框的值，传给search

  // tech-bin
  // const techBin = {
  //   "10x Visium": "55μm",
  //   "Spatial Transcriptomics": "100μm",
  //   "Stereo-Seq": "220nm",
  //   "Slide-seqV2": "10μm",
  //   "Slide-seq": "10μm",
  //   "DBiT-seq": "10μm",
  //   "Seq-Scope": "0.5-0.8μm",
  //   "HDST": "2μm",
  //   "sciMAP-ATAC-seq": "100μm"
  // }

  const techBin = {
    // "2μm": ["Seq-Scope","HDST"],
    // "10μm": ["Slide-seqV2","Slide-seq","DBiT-seq"],
    // "50μm": ["10x Visium"],
    // "100μm": ["Spatial Transcriptomics","sciMAP-ATAC-seq"],
    // "200nm": ["Stereo-Seq"],
    // "50-100μm": [],
    "10μm": ["Seq-Scope","HDST","Stereo-Seq" ],
    "25μm": ["Slide-seqV2","Slide-seq","DBiT-seq"],
    "50μm": ["10x Visium","Spatial Transcriptomics","sciMAP-ATAC-seq"],
    "single-cell SRT": [],
  }

  const cellBinTech = [
    "GeoMx DSP",
    "MERFISH",
    "ClampFISH",
    "EASI-FISH",
    "seqFish",
    "seqFISH+",
    "STARmap",
  ]


  // 有搜索的时候
  const handleSearchComplete = (datas, param) => {
    // datas是从search这个子组件传过来的
    //console.log(datas)
    // data直接改掉src
    setSrc(prevSrc => ({
      ...prevSrc, data: {
        attributes: prevSrc.data.attributes,
        data: datas.data,
        // add lable and item
        label: param.type === 'all' ? 'all' : 'Search',
        item: 'None',
      }
    }))
    TabRef.current?.Fresh(true)
    enterLoading(0, TabRef.current?.Loading)

    // 下面的options都清空
    setSelectedTags(null)
    if (datas.type == 'all') {
      // Rerender('all', 'all')
      // 搜索框清空
      steSearchValue('')
    }
  }

  const handleChange = (tag, checked, label) => {
    // 选择组件时也需要清空搜索框
    steSearchValue('')
    if (checked) {
      setSelectedTags(tag)
      Rerender(label, tag)
    } else {
      setSelectedTags(null)
      Rerender('all', tag)
    }
  }
  //const [selectedKey, setSelectedKey] = useState('datasets')
  //const [visible, setVisible] = useState(true)
  // const onClickTopMenu = (e) => {
  //   setSelectedKey(e.key)
  // }

  const response = useQuery({
    queryKey: ['db'],
    queryFn: () => axios.get('/api/datasets').then((response) => {
      return response.data
    }).catch((error) => {
      console.log(error)
    }),
    staleTime: Infinity,
    retry: false,
    refetchOnWindowFocus: false,
  })

  const filterSpecies = (species) => {
    return (response.data.data.filter((item) => item[3]?.split(";").includes(species)))
  }

  const filterTissues = (tissues) => {
    return (response.data.data.filter((item) => item[4]?.split(";").includes(tissues)))
  }

  const filterTechs = (species, bin) => {
    if(bin){
      const filteredData = response.data.data.filter((item) => {
        const splitArray = item[12]?.split(";");
        if (!splitArray) 
          return false;
        let speciesArray;
        // if(String(species) === "50-100μm"){
        //   let spArray = Technologies.slice(0, 33)
        //   let existTechs = []
        //   Object.entries(techBin).forEach(([key, values]) => {
        //     existTechs.push(values);
        //   });
        //   existTechs = existTechs.flat() // dimension: 2->1
        //   speciesArray = spArray.filter(item => !existTechs.includes(item));
        // }
        if(String(species) == "single-cell SRT"){
          speciesArray = cellBinTech
        }
        else
          speciesArray = techBin[species]
        // console.log(speciesArray)
        return speciesArray.some((species) => splitArray.includes(species));
      })
      return filteredData
    
    }
    else{
      return (response.data.data.filter((item) => item[12]?.split(";").includes(species)))
    }
  }

  const filterDiseases = (disease) => {
    return (response.data.data.filter((item) => item[15]?.toLowerCase().includes(disease.toLowerCase())))
  }

  useEffect(() => {
    if (response.status === 'success' && typeof response.data !== 'undefined') {
      // console.log(response.data)
      setSrc(response)
    }
  }, [response.data])

  const topMenuItems = [
    {
      label: 'Datasets',
      key: 'datasets',
      icon: <DatabaseOutlined />,
    },
    {
      label: 'Samples',
      key: 'samples',
      icon: <ProfileOutlined />,
    },
    {
      label: 'Pairs',
      key: 'pairs',
      icon: <OneToOneOutlined />,
    },
  ]
  const { token } = theme.useToken()
  const TabRef = React.createRef()

  const Rerender = (label, item) => {
    // console.log(label)
    // console.log(item)
    // console.log(src)
    let newSrcData
    if (label === "tissues") {
      newSrcData = filterTissues(item)
      TabRef.current?.GraphConfig(`{tissues: \"${item}\"}`)
    }
    else if (label === "species") {
      newSrcData = filterSpecies(item)
      TabRef.current?.GraphConfig(`{species: \"${item}\"}`)
    }
    else if (label === "technologies") {
      let binFlag = item in techBin? true : false
      newSrcData = filterTechs(item, binFlag)
      TabRef.current?.GraphConfig(`{technologies: \"${item}\"}`)
    }
    else if (label === "disease") {
      newSrcData = filterDiseases(item)
    } else if (label === 'all') {
      newSrcData = response.data.data
    } else {
      newSrcData = src.data.data
    }
    // console.log(newSrcData)
    setSrc(prevSrc => ({
      ...prevSrc, data: {
        attributes: prevSrc.data.attributes,
        data: newSrcData,
        // add lable and item
        label: label,
        item: item,
      }
    }))
    TabRef.current?.Fresh(true)
    enterLoading(0, TabRef.current?.Loading)
  }

  const species = []

  return (
    <Layout
      style={{
        minHeight: '100vh',
      }}>
      {/* <Sider
        width={200}
        breakpoint="lg"
        collapsible
        collapsed={collapsed}
        onCollapse={(value) => setCollapsed(value)}>
        <SideMenu items={srcMeta} collapsed={collapsed} />
      </Sider> */}
      <Layout>
        <Content
          style={{
            padding: '1rem 1rem',
          }}>
          <Col span={6} offset={9} id='Search'>
            <div style={{
              width: '300px',
              height: "70px",
              marginBottom: 30,
              backgroundImage: `url(${logoFig})`,
              objectFit: "cover",
              backgroundSize: 'cover',
              alignItems: "center"
            }} />
          </Col>
          <Col span={16} offset={4}>
            <Search
              value={searchValue}
              onChange={steSearchValue}
              onSearchComplete={handleSearchComplete}
            />
          </Col>
          <Col span={16} offset={4}>
            <Divider style={{ marginTop: 16, marginBottom: 14 }} />
          </Col>
          <Col span={18} offset={3}>
            <TagCollapse
              tags={Species.map((tag) => (<Tag.CheckableTag style={{ margin: 3, fontSize: 16 }}
                key={tag}
                checked={selectedTags === tag}
                onChange={(checked) => {
                  handleChange(tag, checked, 'species')
                }}>{tag}</Tag.CheckableTag>))}
              threshold={5}
              prefix={<b style={{ margin: 3 }}>Species:</b>}
            />
          </Col>
          <Col span={18} offset={3}>
            <TagCollapse
              tags={Tissues.map((tag) => (<Tag.CheckableTag style={{ margin: 3, fontSize: 16 }}
                key={tag}
                checked={selectedTags === tag}
                onChange={(checked) => {
                  handleChange(tag, checked, 'tissues')
                }}>{tag}</Tag.CheckableTag>))}
              threshold={10}
              prefix={<b style={{ margin: 3 }}>Tissues:</b>}
            />
          </Col>
          <Col span={18} offset={3}>
            <TagCollapse
              tags={Technologies.map((tag) => 
                (<Tag.CheckableTag style={{ margin: 3, fontSize: 16 }}
                key={tag}
                checked={selectedTags === tag}
                onChange={(checked) => {
                  handleChange(tag, checked, 'technologies')
                }}>
                  {/* {tag}{techBin[tag]? "("+techBin[tag]+")":''} */}
                  {tag}
                  </Tag.CheckableTag>
                  ))}
              threshold={6}
              prefix={<b style={{ margin: 3 }}>Platforms:</b>}
            />
          </Col>
          <Col span={18} offset={3}>
            <TagCollapse
              tags={Topics.map((tag) => (<Tag.CheckableTag style={{ margin: 3, fontSize: 16 }}
                key={tag}
                checked={selectedTags === tag}
                onChange={(checked) => {
                  handleChange(tag, checked, 'disease')
                }}>{tag}</Tag.CheckableTag >))}
              threshold={6}
              prefix={<b style={{ margin: 3 }}>Topics:</b>}
            />
          </Col>
          <Divider style={{ marginTop: 14, marginBottom: 18 }} />
          <Card>
            <DatasetTab sendData={setSrcMeta} response={src} onRef={TabRef} />
          </Card>
          {/* <div>{JSON.stringify(Array.from(new Set(src.data?.data.map((item) => item[4]))))}</div> */}
        </Content>
        <Footer>
          <Row>
            <Col span={4}>
              <Space direction='horizontal'>
                <h3>Contact us
                  {' '}
                  <img src={contactImg} style={{ width: '20%', height: '20%' }}></img>
                </h3>
              </Space>
              <a href="https://beian.miit.gov.cn" target='_blank'
                style={{
                  color: 'rgba(0, 0, 0, 0.6)',
                  textDecoration: 'none',
                  fontSize: 13
                }}>
                津ICP备2024022781号-1
              </a>
            </Col>
            <Col offset={1}>
              <p>
                *Zhihan Ruan, Centre for Bioinformatics and Intelligent
                Medicine, Nankai University, rrrzhan@nankai.edu.cn
              </p>
              <p>
                *Jian Liu, State Key Laboratory of Medical Chemical Biology, College of Computer Science, Nankai University, jianliu@nankai.edu.cn
              </p>
            </Col>
          </Row>
        </Footer>
      </Layout>
    </Layout>
  )
}
