import Accordion from 'react-bootstrap/Accordion'
import { useAccordionButton } from 'react-bootstrap/AccordionButton'
import Card from 'react-bootstrap/Card'
import { Stack } from 'react-bootstrap'
import PropTypes from 'prop-types'
import { Button } from 'antd'
import { QuestionCircleOutlined } from '@ant-design/icons'

function CustomToggle({ children, eventKey }) {
  const decoratedOnClick = useAccordionButton(eventKey, () =>
    console.log('totally custom!')
  )

  return (
    <QuestionCircleOutlined
      style={{ border: 'secondary', size: 'xxl' }}
      onClick={decoratedOnClick}></QuestionCircleOutlined>
  )
}

const ToggleAccordion = ({ header, children, border }) => {
  return (
    <>
      <Accordion>
        <Card bg="secret" border={border}>
          <Card.Header>
            <Stack direction="horizontal" gap={2}>
              {header}
              <CustomToggle eventKey="0" />
            </Stack>
          </Card.Header>
          <Accordion.Collapse eventKey="0">
            <Card.Body>{children}</Card.Body>
          </Accordion.Collapse>
        </Card>
      </Accordion>
    </>
  )
}

ToggleAccordion.defaultProps = {
  border: null,
}

export default ToggleAccordion
