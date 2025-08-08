import { Component, signal } from '@angular/core';
import { RouterOutlet } from '@angular/router';
import {Map} from './map/map';

@Component({
  selector: 'app-root',
  standalone: true,
  imports: [Map],
  templateUrl: './app.html',
  styleUrl: './app.scss'
})
export class App {
}
