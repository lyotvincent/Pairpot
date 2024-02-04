import React, { useState } from 'react'
import ButtonGroup from 'react-bootstrap/ButtonGroup'
import ToggleButton from 'react-bootstrap/ToggleButton'

export const SwitchButton = () => {
  const [checked, setChecked] = useState(false)
  const [radioValue, setRadioValue] = useState('1')

  const radios = [
    { name: 'dark', value: '1' },
    { name: 'vintage', value: '2' },
    { name: 'macarons', value: '3' },
  ]

  return (
    <>
      <ButtonGroup>
        {radios.map((radio, idx) => (
          <ToggleButton
            key={idx}
            id={`radio-${idx}`}
            type="radio"
            variant="outline-primary"
            name="radio"
            value={radio.value}
            checked={radioValue === radio.value}
            onChange={(e) => setRadioValue(e.currentTarget.value)}>
            {radio.name}
          </ToggleButton>
        ))}
      </ButtonGroup>
    </>
  )
}
