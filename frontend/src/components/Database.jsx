import React, { useState, useEffect, View } from 'react'
import TrackVisibility from 'react-on-screen'
import 'bootstrap/dist/css/bootstrap.min.css'
import { ReactComponent as Female } from '../assets/img/female.svg'
import { ReactComponent as FemaleOrgans } from '../assets/img/female-organs.svg'
import SideMenu from './utils/SideMenu'
import axios from 'axios'
import { useQuery } from 'react-query'
import { Card, Layout, Menu, theme, Space, Button, Tag, Col, Divider, Flex } from 'antd'
import logoFig from "../assets/img/mylogo.png"
import loadingTips from './charts/LoadingTip'
import Loading from './charts/Loading'
import {
  MenuFoldOutlined,
  MenuUnfoldOutlined,
  DatabaseOutlined,
  ProfileOutlined,
  OneToOneOutlined,
} from '@ant-design/icons'
import DatasetTab from './utils/DatasetTab'
import Search from './utils/Search'
import SampleTab from './utils/SampleTab'
const { Header, Content, Footer, Sider } = Layout
const { enterLoading } = Loading

export const Database = () => {
  const [collapsed, setCollapsed] = useState(false)
  const [srcMeta, setSrcMeta] = useState({}) // get the meta data from childrens
  const [src, setSrc] = useState({})
  const [selectedKey, setSelectedKey] = useState('datasets')
  const [visible, setVisible] = useState(true)
  const onClickTopMenu = (e) => {
    setSelectedKey(e.key)
  }

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

  const filterTissues = (species) => {
    return (response.data.data.filter((item) => item[4]?.split(";").includes(species)))
  }

  const filterTechs = (species) => {
    return (response.data.data.filter((item) => item[12]?.split(";").includes(species)))
  }

  const filterDiseases = (diseases) => {
    return (response.data.data.filter((item) => item[16]?.split(" ").includes(diseases)))
  }

  useEffect(() => {
    console.log("toggled")
    if (response.status === 'success' && typeof response.data !== 'undefined') {
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
    let newSrcData
    if(label==="tissues"){
      newSrcData = filterTissues(item)
      TabRef.current?.GraphConfig(`{tissues: \"${item}\"}`)
    }
    else if(label==="species"){
      newSrcData = filterSpecies(item)
    }
    else if(label==="technologies"){
      newSrcData = filterTechs(item)
    }
    else if(label==="diseases"){
      newSrcData = filterDiseases(item)
    } else{
      newSrcData = src.data.data
    }
    setSrc(prevSrc => ({
      ...prevSrc, data: {
        attributes: prevSrc.data.attributes,
        data: newSrcData
      }
    }))
    TabRef.current?.Fresh(true)
    enterLoading(0, TabRef.current?.Loading)
  }

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
          <Col span={16} offset={4}> <Search /> </Col>
          <Col span={16} offset={4}> 
          <Divider orientation="center">Preset Keywords</Divider> 
          </Col>
          <Col span={16} offset={4}> 
          <Flex gap={6} wrap justify='left'>
            Species:
            <Tag color="magenta"
              onClick={() => {
                let newSrcData = filterSpecies("Homo sapiens")
                setSrc(prevSrc => ({
                  ...prevSrc, data: {
                    attributes: prevSrc.data.attributes,
                    data: newSrcData
                  }
                }))
                TabRef.current?.Fresh(true)
                enterLoading(0, TabRef.current?.Loading)
              }}>
              Homo sapiens
            </Tag>
            <Tag color="red" onClick={() => {
              Rerender("species", 'Mus musculus')
            }}>Mus musculus</Tag>
            <Tag color="volcano" onClick={() => {
              Rerender("species", 'Oryctolagus cuniculus')
            }}>
              Oryctolagus cuniculus
            </Tag>
            <Tag color="orange" onClick={() => {
              Rerender("species", 'Danio rerio')
            }}>
              Danio rerio
            </Tag>
            <Tag color="gold" onClick={() => {
              Rerender("species", 'Gallus gallus')
            }}>
              Gallus gallus
            </Tag>
          </Flex>
          </Col>
          <br />
          <Col span={16} offset={4}> 
          <Flex gap={6} wrap justify='left'>
            Tissues:
            <Tag color="magenta"
              onClick={() => {
                Rerender('tissues', 'Brain')
              }}>
              Brain
            </Tag>
            <Tag color="red"
              onClick={() => {
                Rerender('tissues',"Skin")
              }}>
              Skin
            </Tag>
            <Tag color="volcano"
              onClick={() => {
                Rerender('tissues',"Breast")
              }}>
              Breast
            </Tag>
            <Tag color="orange"
              onClick={() => {
                Rerender('tissues',"Kidney")
              }}>
              Kidney
            </Tag>
            <Tag color="gold"
              onClick={() => {
                Rerender('tissues',"Heart")
              }}>
              Heart
            </Tag>
            <Tag color="lime"
              onClick={() => {
                Rerender('tissues',"Lung")
              }}>
              Lung
            </Tag>
            <Tag color="green"
              onClick={() => {
                Rerender('tissues',"Muscle")
              }}>
              Muscle
            </Tag>
            <Tag color="cyan"
              onClick={() => {
                Rerender('tissues',"Prostate")
              }}>
              Prostate
            </Tag>
            <Tag color="blue"
              onClick={() => {
                Rerender('tissues',"Overian")
              }}>
              Overian
            </Tag>
            <Tag color="geekblue"
              onClick={() => {
                Rerender('tissues',"Pancreas")
              }}>
              Pancreas
            </Tag>
            <Tag color="purple"
              onClick={() => {
                Rerender('tissues',"Spinal Cord")
              }}>
              Spinal Cord
            </Tag>
          </Flex>
          </Col>
          <br />
          <Col span={16} offset={4}> 
          <Flex gap={6} wrap justify='left'>
            Technologies:
            <Tag color="magenta"
              onClick={() => {
                Rerender('technologies',"10x Visium")
              }}>
              10x Visium
            </Tag>
            <Tag color="red"
              onClick={() => {
                Rerender('technologies',"scRNA")
              }}>
              scRNA-seq
            </Tag>
            <Tag color="volcano"
              onClick={() => {
                let newSrcData = filterTechs("snRNA")
                setSrc(prevSrc => ({
                  ...prevSrc, data: {
                    attributes: prevSrc.data.attributes,
                    data: newSrcData
                  }
                }))
                TabRef.current?.Fresh(true)
                enterLoading(0, TabRef.current?.Loading)
              }}>
              snRNA-seq
            </Tag>
            <Tag color="orange"
              onClick={() => {
                Rerender('technologies',"Stereo-Seq")
              }}>
              Stereo-Seq
            </Tag>
            <Tag color="gold"
              onClick={() => {
                Rerender('technologies',"MERFISH")
              }}>
              MERFISH
            </Tag>
            <Tag color="lime"
              onClick={() => {
                Rerender('technologies',"Spatial Transcriptomics")
              }}>
              Spatial Transcriptomics
            </Tag>
          </Flex>
          </Col>
          <br/>
          <Col span={16} offset={4}> 
          <Flex gap={6} wrap justify='left'>
            Diseases:
            <Tag color="magenta"
              onClick={() => {
                Rerender("diseases","cancer")
              }}>
              Cancer
            </Tag>
            <Tag color="red"
              onClick={() => {
                Rerender("diseases","injury")
              }}>
              Injury
            </Tag>
            <Tag color="volcano"
              onClick={() => {
                Rerender("diseases","Alzheimer")
              }}>
              Alzheimer
            </Tag>
            <Tag color="orange"
              onClick={() => {
                Rerender("diseases","melanoma")
              }}>
              Melanoma
            </Tag>
            <Tag color="gold"
              onClick={() => {
                Rerender("diseases","carcinoma")
              }}>
              Carcinoma
            </Tag>
          </Flex>
          </Col>
          <Divider />
          <Card>
            <DatasetTab sendData={setSrcMeta} response={src} onRef={TabRef} />
          </Card>
          {/* <div>{JSON.stringify(src)}</div> */}
        </Content>
        <Footer
          style={{
            textAlign: 'center',
          }}>
          PairPot Copyright © 2024. Centre for Bioinformatics and Intelligent Medicine, Nankai University.
        </Footer>
      </Layout>
    </Layout>
  )
}
