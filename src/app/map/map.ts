import {AfterViewInit, ChangeDetectionStrategy, ChangeDetectorRef, Component, inject, signal} from '@angular/core';
import {HttpClient, HttpClientModule} from '@angular/common/http';
import * as L from 'leaflet';
import proj4 from 'proj4';
import {featureToH3Set} from 'geojson2h3';
import * as h3 from 'h3-js';
import {NgIf} from '@angular/common';

@Component({
  selector: 'app-map',
  imports: [HttpClientModule, NgIf],
  standalone: true,
  templateUrl: './map.html',
  styleUrl: './map.scss',
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class Map implements AfterViewInit {
  cdr = inject(ChangeDetectorRef);
  private map!: L.Map;
  private http = inject(HttpClient);
  private geojsonData: any = null;
  private drawnHexagons: L.LayerGroup = L.layerGroup();
  resolution = signal<number>(0);
  zoom = signal<number>(0);
  redrawDebounceTimer: any;
  isLoading = signal<boolean>(true);

  ngAfterViewInit(): void {
    this.initMap();
    this.loadData();
  }

  private initMap(): void {
    this.map = L.map('map').setView([24.7136, 46.6753], 5);

    L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
      maxZoom: 18,
      attribution: '&copy; OpenStreetMap contributors',
    }).addTo(this.map);
    this.zoom.set(this.map.getZoom());

    const scheduleRedraw = () => {
      this.zoom.set(this.map.getZoom());
      this.drawnHexagons.clearLayers();
      if (this.redrawDebounceTimer) clearTimeout(this.redrawDebounceTimer);
      this.isLoading.set(true);

      this.redrawDebounceTimer = setTimeout(() => {
        this.redrawHexagons();
      }, 1000);
    };

    this.map.on('moveend', scheduleRedraw);
    this.map.on('zoomend', scheduleRedraw);

  }

  private loadData(): void {
    this.http.get<any>('assets/data.json').subscribe((geojson) => {
      this.geojsonData = geojson;
      this.redrawHexagons();
    });
  }

  private redrawHexagons(): void {
    if (!this.geojsonData) return;
    this.drawnHexagons.clearLayers();

    const hexLayer = this.drawHexagons(this.geojsonData);
    this.drawnHexagons = hexLayer;
    this.drawnHexagons.addTo(this.map);
  }

  private drawHexagons(geojson: any): L.LayerGroup {
    proj4.defs('EPSG:3857', '+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs');

    this.resolution.set(this.getResolutionFromZoom(this.map.getZoom()));
    this.zoom.set(this.map.getZoom());
    const layerGroup = L.layerGroup();
    const mapBounds = this.map.getBounds();

    geojson.features.forEach((feature: any) => {
      const color = `#${feature.properties.COLOR_HEX}`;

      const converted = JSON.parse(JSON.stringify(feature));
      converted.geometry.coordinates = feature.geometry.coordinates.map((polygon: any[][]) =>
        polygon.map((ring: any[]) =>
          ring.map(([x, y]) => {
            const [lon, lat] = proj4('EPSG:3857', 'EPSG:4326', [x, y]);
            return [lat, lon];
          })
        )
      );


      const h3Indexes: string[] = featureToH3Set(converted, this.resolution());

      h3Indexes.forEach((h3Index) => {
        const hexBoundary = h3.cellToBoundary(h3Index, true);

        const [centerLng, centerLat] = h3.cellToLatLng(h3Index);

        if (mapBounds.contains(L.latLng(centerLat, centerLng))) {
          const latlngs = hexBoundary.map(([lat, lng]) => [lat, lng] as [number, number]);
          const polygon = L.polygon(latlngs, {
            color,
            weight: 1,
            fillOpacity: 0.6,
          });

          layerGroup.addLayer(polygon);
        }
      });
    });
    this.isLoading.set(false);

    return layerGroup;
  }

  private getResolutionFromZoom(zoom: number): number {
    if (zoom <= 5) return 3;
    if (zoom <= 6) return 4;
    if (zoom <= 7) return 5;
    if (zoom <= 8) return 6;
    return 6;
  }
}
